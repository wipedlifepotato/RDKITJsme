"""
Enhanced RDKit API Router — расширенный функционал хемоинформатики
"""

from fastapi import APIRouter, HTTPException, Query, Depends, Response
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from pydantic import BaseModel
from typing import List, Optional
import urllib.parse
import urllib.request

from database import get_db, CachedName
from Molecule import Molecule

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import QED, Draw, inchi
from rdkit.Chem import Lipinski, Crippen, MolSurf, GraphDescriptors
from rdkit.Chem import rdchem

router = APIRouter(prefix="/rdkit/api", tags=["RDKit Enhanced API"])

PROXY_ADDRESS = "127.0.0.1:9050"


@router.get("/", include_in_schema=False)
async def rdkit_index():
    """Индексная страница RDKit API со списком всех эндпоинтов."""
    endpoints = {
        "ADME & Drug-likeness": {
            "GET /rdkit/api/adme": "ADME-свойства, Липинский, QED, SA Score",
        },
        "Analysis": {
            "GET /rdkit/api/functional_groups": "Поиск функциональных групп (25+ типов)",
            "GET /rdkit/api/validate": "Валидация SMILES с диагностикой",
            "GET /rdkit/api/pka": "Оценка pKa (кислотные/основные центры)",
            "GET /rdkit/api/descriptors": "Все дескрипторы RDKit (200+)",
        },
        "Reactions & Conversion": {
            "GET /rdkit/api/reaction": "Применение реакций (SMIRKS)",
            "GET /rdkit/api/by_inchikey": "Поиск по InChIKey",
        },
        "3D & Visualization": {
            "GET /rdkit/api/3d": "3D-структура (SDF/MOL/PDB)",
            "GET /rdkit/api/render_highlighted": "Рендеринг с подсветкой групп",
        },
        "Batch & Library": {
            "POST /rdkit/api/batch": "Пакетный анализ молекул",
            "GET /rdkit/api/library": "Список молекул в библиотеке",
            "POST /rdkit/api/library/save": "Сохранение в библиотеку",
            "DELETE /rdkit/api/library/{smiles}": "Удаление из библиотеки",
        },
        "Comparison": {
            "POST /rdkit/api/compare": "Сравнение нескольких молекул (матрица Tanimoto)",
        },
    }
    return {"endpoints": endpoints, "total": sum(len(v) for v in endpoints.values())}


def set_proxy(proxy: str):
    global PROXY_ADDRESS
    PROXY_ADDRESS = proxy


# ─────────────────────────────────────────────────────────────────────────────
# ADME / Drug-likeness
# ─────────────────────────────────────────────────────────────────────────────


class ADMEResponse(BaseModel):
    smiles: str
    formula: str
    molecular_weight: float
    logp: float
    hbd: int
    hba: int
    tpsa: float
    rotatable_bonds: int
    lipinski_violations: int
    lipinski_pass: bool
    qed: float
    qed_interpretation: str
    sa_score: Optional[float] = None
    num_rings: int
    num_aromatic_rings: int
    num_heteroatoms: int
    num_atoms: int
    num_bonds: int
    fraction_sp3: float
    num_stereocenters: int
    num_undefined_stereocenters: int


@router.get(
    "/adme", response_model=ADMEResponse, summary="ADME / Drug-likeness свойства"
)
def get_adme_properties(smiles: str = Query(..., description="SMILES молекулы")):
    """
    Полный набор ADME-свойств и оценки drug-likeness:
    - Правило Липинского (Rule of Five)
    - QED (Quantitative Estimate of Drug-likeness)
    - SA Score (Synthetic Accessibility)
    - Базовые дескрипторы
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES string")

        mol = m.m

        # Lipinski's Rule of Five
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        violations = 0
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1

        # QED
        qed_score = QED.qed(mol)
        if qed_score >= 0.7:
            qed_interp = "Excellent (drug-like)"
        elif qed_score >= 0.5:
            qed_interp = "Good (moderate drug-likeness)"
        elif qed_score >= 0.3:
            qed_interp = "Fair (poor drug-likeness)"
        else:
            qed_interp = "Poor (unlikely drug-like)"

        # SA Score (synthetic accessibility)
        try:
            from rdkit.Contrib.SA_Score import sascorer

            sa = sascorer.calculateScore(mol)
        except ImportError:
            sa = None

        # Stereocenters
        stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        defined = Chem.FindMolChiralCenters(mol, includeUnassigned=False)

        return ADMEResponse(
            smiles=smiles,
            formula=rdMolDescriptors.CalcMolFormula(mol),
            molecular_weight=round(mw, 4),
            logp=round(logp, 4),
            hbd=hbd,
            hba=hba,
            tpsa=round(Descriptors.TPSA(mol), 4),
            rotatable_bonds=Descriptors.NumRotatableBonds(mol),
            lipinski_violations=violations,
            lipinski_pass=(violations <= 1),
            qed=round(qed_score, 4),
            qed_interpretation=qed_interp,
            sa_score=round(sa, 2) if sa is not None else None,
            num_rings=rdMolDescriptors.CalcNumRings(mol),
            num_aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
            num_heteroatoms=rdMolDescriptors.CalcNumHeteroatoms(mol),
            num_atoms=mol.GetNumAtoms(),
            num_bonds=mol.GetNumBonds(),
            fraction_sp3=round(Lipinski.FractionCSP3(mol), 4),
            num_stereocenters=len(stereo),
            num_undefined_stereocenters=len(stereo) - len(defined),
        )
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Функциональные группы
# ─────────────────────────────────────────────────────────────────────────────

FUNCTIONAL_GROUPS = {
    "hydroxyl": {"smarts": "[OX2H]", "name": "Hydroxyl (-OH)"},
    "carboxylic_acid": {
        "smarts": "[CX3](=O)[OX2H1]",
        "name": "Carboxylic acid (-COOH)",
    },
    "amine_primary": {"smarts": "[NX3;H2;!$(NC=O)]", "name": "Primary amine (-NH2)"},
    "amine_secondary": {
        "smarts": "[NX3;H1;!$(NC=O)]",
        "name": "Secondary amine (-NH-)",
    },
    "amine_tertiary": {"smarts": "[NX3;H0;!$(NC=O)]", "name": "Tertiary amine (-N<)"},
    "amide": {"smarts": "[NX3][CX3](=[OX1])", "name": "Amide (-C(=O)N-)"},
    "nitro": {"smarts": "[NX3](=[OX1])(=[OX1])", "name": "Nitro (-NO2)"},
    "sulfonyl": {"smarts": "[SX4](=[OX1])(=[OX1])", "name": "Sulfonyl (-SO2-)"},
    "phosphate": {"smarts": "[PX4](=[OX1])([OX2])[OX2]", "name": "Phosphate"},
    "aldehyde": {"smarts": "[CX3H1](=O)[#6]", "name": "Aldehyde (-CHO)"},
    "ketone": {"smarts": "[#6][CX3](=O)[#6]", "name": "Ketone (>C=O)"},
    "ester": {"smarts": "[#6][CX3](=O)[OX2H0][#6]", "name": "Ester (-COOR)"},
    "ether": {"smarts": "[OD2]([#6])[#6]", "name": "Ether (-O-)"},
    "alkene": {"smarts": "[CX3]=[CX3]", "name": "Alkene (C=C)"},
    "alkyne": {"smarts": "[CX2]#[CX2]", "name": "Alkyne (C≡C)"},
    "aromatic_ring": {"smarts": "c1ccccc1", "name": "Aromatic ring (benzene)"},
    "heterocycle": {"smarts": "[!#6]1~[!#6]~*~*~*~1", "name": "Heterocycle"},
    "halogen_f": {"smarts": "[F]", "name": "Fluorine"},
    "halogen_cl": {"smarts": "[Cl]", "name": "Chlorine"},
    "halogen_br": {"smarts": "[Br]", "name": "Bromine"},
    "halogen_i": {"smarts": "[I]", "name": "Iodine"},
    "thiol": {"smarts": "[SX2H]", "name": "Thiol (-SH)"},
    "nitrile": {"smarts": "[CX2]#[NX1]", "name": "Nitrile (-C≡N)"},
    "azide": {"smarts": "[NX2]=[NX2]=[NX1]", "name": "Azide (-N3)"},
    "imine": {"smarts": "[CX3]=[NX2]", "name": "Imine (C=N)"},
    "hydrazine": {"smarts": "[NX3][NX3]", "name": "Hydrazine (-N-N-)"},
    "urea": {"smarts": "[NX3][CX3](=[OX1])[NX3]", "name": "Urea"},
    "carbamate": {"smarts": "[OX2][CX3](=[OX1])[NX3]", "name": "Carbamate"},
    "sulfonamide": {"smarts": "[SX4](=[OX1])(=[OX1])[NX3]", "name": "Sulfonamide"},
}


@router.get("/functional_groups", summary="Поиск функциональных групп")
def find_functional_groups(smiles: str = Query(..., description="SMILES молекулы")):
    """
    Определяет наличие функциональных групп в молекуле по SMARTS-паттернам.
    Возвращает список найденных групп с количеством совпадений.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES string")

        mol = m.m
        found = []
        for key, info in FUNCTIONAL_GROUPS.items():
            pattern = Chem.MolFromSmarts(info["smarts"])
            if pattern is None:
                continue
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                found.append(
                    {
                        "key": key,
                        "name": info["name"],
                        "smarts": info["smarts"],
                        "count": len(matches),
                        "atom_indices": [list(m) for m in matches],
                    }
                )

        return {
            "smiles": smiles,
            "total_groups_found": len(found),
            "functional_groups": found,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Валидация SMILES с диагностикой
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/validate", summary="Валидация SMILES с подробной диагностикой")
def validate_smiles(smiles: str = Query(..., description="SMILES для проверки")):
    """
    Проверяет корректность SMILES и возвращает подробную диагностику:
    - Тип ошибки (valence, aromaticity, syntax)
    - Позицию ошибки
    - Предложения по исправлению
    """
    result = {
        "smiles": smiles,
        "valid": False,
        "canonical_smiles": None,
        "error_type": None,
        "error_message": None,
        "error_position": None,
        "warnings": [],
    }

    if not smiles or not smiles.strip():
        result["error_type"] = "empty"
        result["error_message"] = "SMILES string is empty"
        return result

    smiles = smiles.strip()

    # Проверка на недопустимые символы
    allowed_chars = set(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]=#@+-./\\:%"
    )
    invalid_chars = set()
    for i, ch in enumerate(smiles):
        if ch not in allowed_chars:
            invalid_chars.add((i, ch))
    if invalid_chars:
        result["error_type"] = "invalid_characters"
        result["error_message"] = f"Invalid characters found: {invalid_chars}"
        result["error_position"] = list(invalid_chars)[0][0]
        return result

    # Проверка баланса скобок
    open_paren = smiles.count("(")
    close_paren = smiles.count(")")
    if open_paren != close_paren:
        result["error_type"] = "unbalanced_parentheses"
        result["error_message"] = (
            f"Unbalanced parentheses: {open_paren} open vs {close_paren} close"
        )
        return result

    open_bracket = smiles.count("[")
    close_bracket = smiles.count("]")
    if open_bracket != close_bracket:
        result["error_type"] = "unbalanced_brackets"
        result["error_message"] = (
            f"Unbalanced brackets: {open_bracket} open vs {close_bracket} close"
        )
        return result

    # Попытка парсинга
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # Попытка получить более детальную ошибку
        try:
            from rdkit.Chem import rdMolDescriptors

            # Пробуем с sanitize=False
            mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol_raw is not None:
                try:
                    Chem.SanitizeMol(mol_raw)
                except Exception as sanitize_err:
                    err_str = str(sanitize_err)
                    if "valence" in err_str.lower():
                        result["error_type"] = "valence_error"
                    elif "aromatic" in err_str.lower():
                        result["error_type"] = "aromaticity_error"
                    elif "kekul" in err_str.lower():
                        result["error_type"] = "kekulization_error"
                    else:
                        result["error_type"] = "sanitization_error"
                    result["error_message"] = err_str
            else:
                result["error_type"] = "syntax_error"
                result["error_message"] = "RDKit could not parse this SMILES string"
        except Exception:
            result["error_type"] = "syntax_error"
            result["error_message"] = "RDKit could not parse this SMILES string"
        return result

    # Валидный SMILES
    result["valid"] = True
    result["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)

    # Предупреждения
    if mol.GetNumAtoms() > 100:
        result["warnings"].append("Large molecule (>100 atoms)")
    if Descriptors.MolWt(mol) > 1000:
        result["warnings"].append("High molecular weight (>1000 Da)")
    if Descriptors.MolLogP(mol) > 7:
        result["warnings"].append("Very high lipophilicity (LogP > 7)")

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Реакции (SMIRKS)
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/reaction", summary="Применение реакции (SMIRKS)")
def apply_reaction(
    smiles: str = Query(..., description="SMILES реагента"),
    reaction_smarts: str = Query(..., description="SMARTS реакции (SMIRKS)"),
):
    """
    Применяет реакцию (SMIRKS) к молекуле.
    Возвращает список продуктов реакции.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")

        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            raise ValueError("Invalid SMIRKS reaction pattern")

        products = rxn.RunReactants((mol,))

        if not products:
            return {
                "smiles": smiles,
                "reaction_smarts": reaction_smarts,
                "products": [],
                "product_count": 0,
                "message": "No reaction products formed",
            }

        unique_products = []
        seen = set()
        for product_set in products:
            for product in product_set:
                try:
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product, canonical=True)
                    if smi not in seen:
                        seen.add(smi)
                        unique_products.append(smi)
                except Exception:
                    continue

        return {
            "smiles": smiles,
            "reaction_smarts": reaction_smarts,
            "products": unique_products,
            "product_count": len(unique_products),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# pKa (оценка)
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/pka", summary="Оценка pKa молекулы")
def estimate_pka(smiles: str = Query(..., description="SMILES молекулы")):
    """
    Оценивает кислотные/основные центры молекулы.
    Использует эвристический подход на основе функциональных групп.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES string")

        mol = m.m
        pka_sites = []

        # Карбоновые кислоты (pKa ~4-5)
        acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        if acid_pattern:
            for match in mol.GetSubstructMatches(acid_pattern):
                pka_sites.append(
                    {
                        "atom_idx": match[0],
                        "type": "carboxylic_acid",
                        "estimated_pka": 4.5,
                        "strength": "acidic",
                    }
                )

        # Фенолы (pKa ~9-10)
        phenol_pattern = Chem.MolFromSmarts("[OX2H]c1ccccc1")
        if phenol_pattern:
            for match in mol.GetSubstructMatches(phenol_pattern):
                pka_sites.append(
                    {
                        "atom_idx": match[0],
                        "type": "phenol",
                        "estimated_pka": 10.0,
                        "strength": "acidic",
                    }
                )

        # Алифатические амины (pKa ~9-11)
        amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O);!$([NH2]c1ccccc1)]")
        if amine_pattern:
            for match in mol.GetSubstructMatches(amine_pattern):
                pka_sites.append(
                    {
                        "atom_idx": match[0],
                        "type": "aliphatic_amine",
                        "estimated_pka": 10.0,
                        "strength": "basic",
                    }
                )

        # Анилины (pKa ~4-5)
        aniline_pattern = Chem.MolFromSmarts("[NX3;H2,H1]c1ccccc1")
        if aniline_pattern:
            for match in mol.GetSubstructMatches(aniline_pattern):
                pka_sites.append(
                    {
                        "atom_idx": match[0],
                        "type": "aniline",
                        "estimated_pka": 4.6,
                        "strength": "basic",
                    }
                )

        # Тиолы (pKa ~8-10)
        thiol_pattern = Chem.MolFromSmarts("[SX2H]")
        if thiol_pattern:
            for match in mol.GetSubstructMatches(thiol_pattern):
                pka_sites.append(
                    {
                        "atom_idx": match[0],
                        "type": "thiol",
                        "estimated_pka": 9.0,
                        "strength": "acidic",
                    }
                )

        # Сульфоновые кислоты (pKa ~-2)
        sulfonic_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])[OX2H1]")
        if sulfonic_pattern:
            for match in mol.GetSubstructMatches(sulfonic_pattern):
                pka_sites.append(
                    {
                        "atom_idx": match[0],
                        "type": "sulfonic_acid",
                        "estimated_pka": -2.0,
                        "strength": "strong_acid",
                    }
                )

        return {
            "smiles": smiles,
            "pka_sites": pka_sites,
            "total_sites": len(pka_sites),
            "acidic_sites": len(
                [s for s in pka_sites if s["strength"] in ("acidic", "strong_acid")]
            ),
            "basic_sites": len([s for s in pka_sites if s["strength"] == "basic"]),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Batch API
# ─────────────────────────────────────────────────────────────────────────────


class BatchRequest(BaseModel):
    smiles_list: List[str]
    properties: Optional[List[str]] = None  # Если None — все свойства


@router.post("/batch", summary="Пакетный анализ молекул")
def batch_analyze(request: BatchRequest):
    """
    Анализирует список SMILES за один запрос.
    Возвращает свойства для каждой молекулы.
    """
    results = []
    for smiles in request.smiles_list:
        try:
            m = Molecule(smiles)
            if m.m is None:
                results.append(
                    {"smiles": smiles, "valid": False, "error": "Invalid SMILES"}
                )
                continue

            mol = m.m
            props = {
                "smiles": smiles,
                "valid": True,
                "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "molecular_weight": round(Descriptors.MolWt(mol), 4),
                "logp": round(Descriptors.MolLogP(mol), 4),
                "hbd": Descriptors.NumHDonors(mol),
                "hba": Descriptors.NumHAcceptors(mol),
                "tpsa": round(Descriptors.TPSA(mol), 4),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "qed": round(QED.qed(mol), 4),
                "num_atoms": mol.GetNumAtoms(),
                "num_bonds": mol.GetNumBonds(),
                "num_rings": rdMolDescriptors.CalcNumRings(mol),
            }

            # Lipinski
            violations = 0
            if props["molecular_weight"] > 500:
                violations += 1
            if props["logp"] > 5:
                violations += 1
            if props["hbd"] > 5:
                violations += 1
            if props["hba"] > 10:
                violations += 1
            props["lipinski_violations"] = violations
            props["lipinski_pass"] = violations <= 1

            results.append(props)
        except Exception as e:
            results.append({"smiles": smiles, "valid": False, "error": str(e)})

    return {
        "total": len(request.smiles_list),
        "valid_count": sum(1 for r in results if r.get("valid")),
        "invalid_count": sum(1 for r in results if not r.get("valid")),
        "results": results,
    }


# ─────────────────────────────────────────────────────────────────────────────
# InChIKey поддержка
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/by_inchikey", summary="Получить свойства по InChIKey")
def get_by_inchikey(inchikey: str = Query(..., description="InChIKey молекулы")):
    """
    Получает свойства молекулы по InChIKey.
    Использует Cactus/NCI resolver для конвертации.
    """
    try:
        # Проверка формата InChIKey
        if len(inchikey) != 27 or inchikey.count("-") != 2:
            raise ValueError(
                "Invalid InChIKey format. Expected format: XXXXXXXXXXXXXX-XXXXXXXXXX-X"
            )

        # Попытка разрешить через Cactus
        safe_key = urllib.parse.quote(inchikey, safe="")
        target_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_key}/smiles"

        opener = urllib.request.build_opener()
        if PROXY_ADDRESS.strip():
            proxy_host, proxy_port = PROXY_ADDRESS.split(":")
            import socks
            from sockshandler import SocksiPyHandler

            proxy_handler = SocksiPyHandler(socks.SOCKS5, proxy_host, int(proxy_port))
            opener = urllib.request.build_opener(proxy_handler)

        req = urllib.request.Request(target_url, headers={"User-Agent": "Mozilla/5.0"})
        with opener.open(req, timeout=10) as resp:
            smiles = resp.read().decode("utf-8").strip()

        if not smiles:
            raise ValueError("Could not resolve InChIKey")

        # Получаем свойства
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Resolved SMILES is invalid")

        mol = m.m
        return {
            "inchikey": inchikey,
            "smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(Descriptors.MolWt(mol), 4),
            "logp": round(Descriptors.MolLogP(mol), 4),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": round(Descriptors.TPSA(mol), 4),
            "qed": round(QED.qed(mol), 4),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Полный набор дескрипторов
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/descriptors", summary="Все дескрипторы RDKit")
def get_all_descriptors(smiles: str = Query(..., description="SMILES молекулы")):
    """
    Возвращает все доступные дескрипторы RDKit для молекулы.
    Включает 200+ молекулярных дескрипторов.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES string")

        mol = m.m

        # Вычисляем все дескрипторы
        descriptors = {}
        for name, func in Descriptors.descList:
            try:
                value = func(mol)
                if isinstance(value, (int, float)):
                    descriptors[name] = (
                        round(value, 6) if isinstance(value, float) else value
                    )
            except Exception:
                descriptors[name] = None

        # Дополнительные дескрипторы
        extra = {
            "NumAtoms": mol.GetNumAtoms(),
            "NumBonds": mol.GetNumBonds(),
            "NumHeavyAtoms": mol.GetNumHeavyAtoms(),
            "NumRings": rdMolDescriptors.CalcNumRings(mol),
            "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "NumAliphaticRings": rdMolDescriptors.CalcNumAliphaticRings(mol),
            "NumSaturatedRings": rdMolDescriptors.CalcNumSaturatedRings(mol),
            "NumHeteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
            "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
            "NumHDonors": Descriptors.NumHDonors(mol),
            "NumHAcceptors": Descriptors.NumHAcceptors(mol),
            "FractionCSP3": Lipinski.FractionCSP3(mol),
            "MolWt": Descriptors.MolWt(mol),
            "ExactMolWt": Descriptors.ExactMolWt(mol),
            "MolLogP": Descriptors.MolLogP(mol),
            "MolMR": Descriptors.MolMR(mol),
            "TPSA": Descriptors.TPSA(mol),
            "BalabanJ": GraphDescriptors.BalabanJ(mol),
            "BertzCT": GraphDescriptors.BertzCT(mol),
            "Chi0": GraphDescriptors.Chi0(mol),
            "Chi1": GraphDescriptors.Chi1(mol),
            "Kappa1": GraphDescriptors.Kappa1(mol),
            "Kappa2": GraphDescriptors.Kappa2(mol),
            "Kappa3": GraphDescriptors.Kappa3(mol),
            "HallKierAlpha": Descriptors.HallKierAlpha(mol),
            "NumValenceElectrons": Descriptors.NumValenceElectrons(mol),
            "NumRadicalElectrons": Descriptors.NumRadicalElectrons(mol),
            "MaxPartialCharge": Descriptors.MaxPartialCharge(mol),
            "MinPartialCharge": Descriptors.MinPartialCharge(mol),
            "MaxAbsPartialCharge": Descriptors.MaxAbsPartialCharge(mol),
            "MinAbsPartialCharge": Descriptors.MinAbsPartialCharge(mol),
        }

        # QED components
        try:
            qed_props = QED.properties(mol)
            extra["QED"] = {
                "MW": round(qed_props.MW, 4),
                "ALOGP": round(qed_props.ALOGP, 4),
                "HBA": qed_props.HBA,
                "HBD": qed_props.HBD,
                "PSA": round(qed_props.PSA, 4),
                "ROTB": qed_props.ROTB,
                "AROM": qed_props.AROM,
                "ALERTS": qed_props.ALERTS,
                "score": round(QED.qed(mol), 4),
            }
        except Exception:
            pass

        return {
            "smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
            "descriptor_count": len(descriptors) + len(extra),
            "descriptors": descriptors,
            "extra_descriptors": extra,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# 3D визуализация (base64 SDF + координаты)
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/3d", summary="3D-структура молекулы")
def get_3d_structure(
    smiles: str = Query(..., description="SMILES молекулы"),
    format: str = Query("sdf", description="Формат: sdf, mol, pdb"),
):
    """
    Генерирует 3D-структуру молекулы.
    Возвращает в формате SDF, MOL или PDB.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = Chem.AddHs(m.m)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)

        if format.lower() == "sdf":
            block = Chem.MolToMolBlock(mol)
            media_type = "chemical/x-mdl-sdfile"
        elif format.lower() == "mol":
            block = Chem.MolToMolBlock(mol)
            media_type = "chemical/x-mdl-molfile"
        elif format.lower() == "pdb":
            block = Chem.MolToPDBBlock(mol)
            media_type = "chemical/x-pdb"
        else:
            raise ValueError(f"Unsupported format: {format}")

        return Response(content=block, media_type=media_type)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# История / Библиотека молекул
# ─────────────────────────────────────────────────────────────────────────────


class MoleculeEntry(BaseModel):
    smiles: str
    name: Optional[str] = None
    tags: Optional[List[str]] = None


@router.post("/library/save", summary="Сохранить молекулу в библиотеку")
def save_to_library(entry: MoleculeEntry, db: Session = Depends(get_db)):
    """Сохраняет молекулу в локальную библиотеку (SQLite)."""
    try:
        m = Molecule(entry.smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        canonical = Chem.MolToSmiles(m.m, canonical=True)

        # Проверяем, есть ли уже
        existing = db.query(CachedName).filter(CachedName.smiles == canonical).first()
        if existing:
            if entry.name:
                existing.name = entry.name
            db.commit()
            return {"status": "updated", "smiles": canonical}

        new_entry = CachedName(smiles=canonical, name=entry.name or "")
        db.add(new_entry)
        db.commit()
        return {"status": "saved", "smiles": canonical}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/library", summary="Список молекул в библиотеке")
def list_library(
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
):
    """Возвращает список сохранённых молекул."""
    try:
        entries = db.query(CachedName).offset(offset).limit(limit).all()
        total = db.query(CachedName).count()
        return {
            "total": total,
            "limit": limit,
            "offset": offset,
            "entries": [{"smiles": e.smiles, "name": e.name} for e in entries],
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.delete("/library/{smiles}", summary="Удалить молекулу из библиотеки")
def delete_from_library(smiles: str, db: Session = Depends(get_db)):
    """Удаляет молекулу из библиотеки по SMILES."""
    try:
        entry = db.query(CachedName).filter(CachedName.smiles == smiles).first()
        if not entry:
            raise HTTPException(status_code=404, detail="Molecule not found in library")
        db.delete(entry)
        db.commit()
        return {"status": "deleted", "smiles": smiles}
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Сравнение нескольких молекул
# ─────────────────────────────────────────────────────────────────────────────


@router.post("/compare", summary="Сравнение нескольких молекул")
def compare_molecules(
    smiles_list: List[str] = Query(..., description="Список SMILES для сравнения")
):
    """
    Сравнивает несколько молекул попарно.
    Возвращает матрицу сходства Tanimoto.
    """
    try:
        mols = []
        valid_smiles = []
        for smiles in smiles_list:
            m = Molecule(smiles)
            if m.m is not None:
                mols.append(m.m)
                valid_smiles.append(smiles)

        if len(mols) < 2:
            raise ValueError("Need at least 2 valid SMILES for comparison")

        fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]

        n = len(mols)
        matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                matrix[i][j] = round(DataStructs.TanimotoSimilarity(fps[i], fps[j]), 4)

        return {"smiles": valid_smiles, "count": n, "similarity_matrix": matrix}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Рендеринг с подсветкой функциональных групп
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/render_highlighted", summary="Рендеринг с подсветкой групп")
def render_highlighted(
    smiles: str = Query(..., description="SMILES молекулы"),
    highlight_groups: Optional[List[str]] = Query(
        None, description="Ключи групп для подсветки"
    ),
):
    """
    Рендерит молекулу с подсветкой функциональных групп.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = m.m
        highlight_atoms = []
        highlight_colors = {}

        if highlight_groups:
            for group_key in highlight_groups:
                if group_key in FUNCTIONAL_GROUPS:
                    pattern = Chem.MolFromSmarts(FUNCTIONAL_GROUPS[group_key]["smarts"])
                    if pattern:
                        matches = mol.GetSubstructMatches(pattern)
                        for match in matches:
                            for atom_idx in match:
                                if atom_idx not in highlight_atoms:
                                    highlight_atoms.append(atom_idx)
                                    highlight_colors[atom_idx] = (
                                        0.2,
                                        0.8,
                                        0.2,
                                    )  # Зелёный

        img = Draw.MolToImage(
            mol,
            size=(500, 500),
            highlightAtoms=highlight_atoms if highlight_atoms else None,
            highlightAtomColors=highlight_colors if highlight_colors else None,
        )

        from io import BytesIO
        import base64

        buffer = BytesIO()
        img.save(buffer, format="PNG")
        b64_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

        return {
            "smiles": smiles,
            "image_base64": b64_str,
            "highlighted_atoms": highlight_atoms,
            "highlighted_groups": highlight_groups or [],
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# ЯМР (1H и 13C) — эвристическая оценка
# ─────────────────────────────────────────────────────────────────────────────


# Таблицы химических сдвигов для 1H ЯМР (ppm)
NMR_1H_SHIFTS = {
    "methyl": {"range": (0.8, 1.2), "desc": "CH3 (метил)"},
    "methylene": {"range": (1.2, 1.5), "desc": "CH2 (метилен)"},
    "methine": {"range": (1.4, 1.7), "desc": "CH (метин)"},
    "allylic": {"range": (1.6, 2.6), "desc": "Аллильный CH"},
    "alpha_to_carbonyl": {"range": (2.0, 2.5), "desc": "CH рядом с C=O"},
    "alkyne": {"range": (2.0, 3.0), "desc": "Ацетиленовый CH"},
    "methoxy": {"range": (3.3, 4.0), "desc": "OCH3 (метокси)"},
    "alpha_to_oxygen": {"range": (3.3, 4.5), "desc": "CH рядом с O"},
    "alkene": {"range": (4.5, 6.5), "desc": "Алкеновый CH"},
    "aromatic": {"range": (6.5, 8.5), "desc": "Ароматический CH"},
    "aldehyde": {"range": (9.0, 10.0), "desc": "Альдегидный CH"},
    "carboxylic_acid": {"range": (10.0, 13.0), "desc": "COOH (карбоновая кислота)"},
    "phenol": {"range": (4.5, 7.0), "desc": "OH (фенол)"},
    "amine": {"range": (1.0, 5.0), "desc": "NH (амин)"},
}

# Таблицы химических сдвигов для 13C ЯМР (ppm)
NMR_13C_SHIFTS = {
    "methyl": {"range": (5, 30), "desc": "CH3 (метил)"},
    "methylene": {"range": (15, 45), "desc": "CH2 (метилен)"},
    "methine": {"range": (25, 50), "desc": "CH (метин)"},
    "alkyne": {"range": (65, 90), "desc": "C≡C (ацетилен)"},
    "alkene": {"range": (100, 150), "desc": "C=C (алкен)"},
    "aromatic": {"range": (110, 160), "desc": "Ароматический C"},
    "carbonyl_ester": {"range": (160, 175), "desc": "C=O (сложный эфир)"},
    "carbonyl_amide": {"range": (160, 180), "desc": "C=O (амид)"},
    "carbonyl_ketone": {"range": (190, 220), "desc": "C=O (кетон, альдегид)"},
    "nitrile": {"range": (115, 125), "desc": "C≡N (нитрил)"},
}


# ─────────────────────────────────────────────────────────────────────────────
# ЯМР (1H и 13C) — улучшенная эвристическая оценка
# ─────────────────────────────────────────────────────────────────────────────

@router.get("/nmr", summary="Оценка химических сдвигов ЯМР (1H и 13C)")
def estimate_nmr(smiles: str = Query(..., description="SMILES молекулы")):
    """
    Эвристическая оценка химических сдвигов ЯМР с учетом химического окружения.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = Chem.AddHs(m.m)

        # 1H ЯМР: анализ по функциональным центрам и метильным/метиленовым группам
        proton_shifts = []
        
        # 1. Поиск карбоновых кислот (-COOH)
        acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        if acid_pattern:
            for match in mol.GetSubstructMatches(acid_pattern):
                h_idx = match[2] # Индекс протона у кислорода в базовом патерне или найдем через соседи
                oxy_atom = mol.GetAtomWithIdx(match[2])
                for n in oxy_atom.GetNeighbors():
                    if n.GetAtomicNum() == 1:
                        proton_shifts.append({
                            "group": "COOH",
                            "shift_range": [10.5, 12.5],
                            "description": "COOH (карбоновая кислота)",
                            "multiplicity": "s",
                            "integration": 1
                        })

        # 2. Поиск метильных групп (-CH3)
        methyl_pattern = Chem.MolFromSmarts("[CX4H3]")
        if methyl_pattern:
            for match in mol.GetSubstructMatches(methyl_pattern):
                c_atom = mol.GetAtomWithIdx(match[0])
                # Проверим, примыкает ли к карбонилу (например, ацетаты/уксусная кислота)
                is_near_carbonyl = False
                for n in c_atom.GetNeighbors():
                    if n.GetAtomicNum() == 6:
                        for nn in n.GetNeighbors():
                            if nn.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n.GetIdx(), nn.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                is_near_carbonyl = True
                
                if is_near_carbonyl:
                    shift_range = [2.0, 2.3]
                    desc = "CH3 соседний с C=O (ацетил/уксусная кислота)"
                else:
                    shift_range = [0.8, 1.2]
                    desc = "CH3 (алифатический метил)"

                # Собираем сами атомы водорода для этого метила
                h_count = sum(1 for n in c_atom.GetNeighbors() if n.GetAtomicNum() == 1)
                if h_count > 0:
                    proton_shifts.append({
                        "group": "CH3",
                        "shift_range": shift_range,
                        "description": desc,
                        "multiplicity": "s" if is_near_carbonyl else "t",
                        "integration": h_count
                    })

        # 13C ЯМР
        carbon_shifts = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Углерод
                c_shift = _estimate_13c_shift_improved(mol, atom)
                if c_shift:
                    carbon_shifts.append({
                        "atom_idx": atom.GetIdx(),
                        "shift_range": c_shift["range"],
                        "description": c_shift["desc"],
                    })

        return {
            "smiles": smiles,
            "proton_nmr": {
                "shifts": proton_shifts,
            },
            "carbon_nmr": {
                "shifts": carbon_shifts,
            },
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


def _estimate_13c_shift_improved(mol, c_atom):
    """Точная оценка сдвига для 13C ЯМР."""
    c_idx = c_atom.GetIdx()
    
    # Проверяем карбонильные группы
    for neighbor in c_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
            # Проверяем, кислота ли это / сложный эфир / кетон
            is_acid_or_ester = False
            for n2 in c_atom.GetNeighbors():
                if n2.GetAtomicNum() == 8 and n2.GetIdx() != neighbor.GetIdx():
                    is_acid_or_ester = True
            if is_acid_or_ester:
                return {"range": [165, 185], "desc": "C=O (карбоновая кислота / сложный эфир)"}
            return {"range": [190, 220], "desc": "C=O (кетон / альдегид)"}

    # Ароматические углероды
    if c_atom.GetIsAromatic():
        return {"range": [110, 160], "desc": "Ароматический C"}

    # Алифатические углероды по количеству водородов
    num_h = sum(1 for n in c_atom.GetNeighbors() if n.GetAtomicNum() == 1)
    if num_h >= 3:
        return {"range": [10, 30], "desc": "CH3 (метил)"}
    elif num_h == 2:
        return {"range": [15, 45], "desc": "CH2 (метилен)"}
    else:
        return {"range": [25, 60], "desc": "CH (метин) или четвертичный углерод"}

def _estimate_1h_shift(mol, h_atom, heavy_atom):
    """Оценка химического сдвига протона."""
    heavy_sym = heavy_atom.GetSymbol()
    heavy_idx = heavy_atom.GetIdx()

    # Проверяем окружение тяжёлого атома
    if heavy_sym == "C":
        # Проверяем, является ли углерод карбонильным
        for neighbor in heavy_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(
                heavy_idx, neighbor.GetIdx()
            ).GetBondType() == Chem.BondType.DOUBLE:
                return NMR_1H_SHIFTS["alpha_to_carbonyl"]

        # Проверяем ароматичность
        if heavy_atom.GetIsAromatic():
            return NMR_1H_SHIFTS["aromatic"]

        # Проверяем алкены
        for neighbor in heavy_atom.GetNeighbors():
            if mol.GetBondBetweenAtoms(
                heavy_idx, neighbor.GetIdx()
            ).GetBondType() == Chem.BondType.DOUBLE:
                return NMR_1H_SHIFTS["alkene"]

        # Определяем тип углерода по числу соседей
        num_h_neighbors = sum(1 for n in heavy_atom.GetNeighbors() if n.GetAtomicNum() == 1)
        if num_h_neighbors >= 3:
            return NMR_1H_SHIFTS["methyl"]
        elif num_h_neighbors == 2:
            return NMR_1H_SHIFTS["methylene"]
        else:
            return NMR_1H_SHIFTS["methine"]

    elif heavy_sym == "O":
        # Проверяем, является ли кислород частью карбоновой кислоты
        for neighbor in heavy_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                for n2 in neighbor.GetNeighbors():
                    if n2.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(
                        neighbor.GetIdx(), n2.GetIdx()
                    ).GetBondType() == Chem.BondType.DOUBLE:
                        return NMR_1H_SHIFTS["carboxylic_acid"]
        # Проверяем фенол
        if heavy_atom.GetIsAromatic():
            return NMR_1H_SHIFTS["phenol"]
        return NMR_1H_SHIFTS["phenol"]

    elif heavy_sym == "N":
        return NMR_1H_SHIFTS["amine"]

    return None

def _estimate_multiplicity(h_atom, heavy_atom):
    """Оценка мультиплетности (n+1 правило)."""
    # Считаем соседние протоны на соседних атомах
    n_protons = 0
    for neighbor in heavy_atom.GetNeighbors():
        if neighbor.GetIdx() == h_atom.GetIdx():
            continue
        if neighbor.GetAtomicNum() == 1:
            n_protons += 1
        elif neighbor.GetAtomicNum() == 6:
            # Считаем протоны на соседнем углероде
            for n2 in neighbor.GetNeighbors():
                if n2.GetAtomicNum() == 1:
                    n_protons += 1

    if n_protons == 0:
        return "s"  # Синглет
    elif n_protons == 1:
        return "d"  # Дублет
    elif n_protons == 2:
        return "t"  # Триплет
    elif n_protons == 3:
        return "q"  # Квартет
    else:
        return "m"  # Мультиплет


# ─────────────────────────────────────────────────────────────────────────────
# Токсичность (ADMET)
# ─────────────────────────────────────────────────────────────────────────────


# Структуры, связанные с токсичностью
TOXICITY_ALERTS = {
    "hERG": {
        "name": "hERG Channel Blockade",
        "description": "Блокирование калиевых каналов hERG — кардиотоксичность",
        "smarts": [
            {"pattern": "[NX3;H2,H1;!$(NC=O)]", "desc": "Основные амины"},
            {"pattern": "[NX3;H2,H1]c1ccccc1", "desc": "Анилины"},
            {"pattern": "c1ccc2c(c1)ccc2", "desc": "Полициклические ароматические углеводороды"},
        ],
        "risk": "Высокий риск удлинения интервала QT",
    },
    "ames": {
        "name": "Ames Test (Mutagenicity)",
        "description": "Мутагенность — повреждение ДНК",
        "smarts": [
            {"pattern": "[NX2]=[NX2]=[NX1]", "desc": "Азиды"},
            {"pattern": "[NX3](=[OX1])(=[OX1])", "desc": "Нитросоединения"},
            {"pattern": "[CX3](=O)[NX2][CX3](=O)", "desc": "Диамиды"},
            {"pattern": "[SX4](=[OX1])(=[OX1])[OX2H1]", "desc": "Сульфоновые кислоты"},
            {"pattern": "[Cl,Br,I][CX4]", "desc": "Галогеналканы"},
        ],
        "risk": "Потенциальные мутагены",
    },
    "hepatotoxicity": {
        "name": "Hepatotoxicity",
        "description": "Токсичность для печени",
        "smarts": [
            {"pattern": "[NX3;H2]c1ccccc1", "desc": "Анилины"},
            {"pattern": "[CX3](=O)[NX2][CX3](=O)", "desc": "Диамиды"},
            {"pattern": "[SX4](=[OX1])(=[OX1])[NX3]", "desc": "Сульфонамиды"},
            {"pattern": "[PX4](=[OX1])([OX2])[OX2]", "desc": "Фосфаты"},
        ],
        "risk": "Потенциальная гепатотоксичность",
    },
}


@router.get("/toxicity", summary="Оценка токсичности (hERG, Ames, гепатотоксичность)")
def assess_toxicity(smiles: str = Query(..., description="SMILES молекулы")):
    """
    Оценка потенциальной токсичности молекулы:
    - hERG (кардиотоксичность)
    - Ames test (мутагенность)
    - Гепатотоксичность
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = m.m
        alerts = {}

        for tox_type, info in TOXICITY_ALERTS.items():
            matches = []
            for smarts_info in info["smarts"]:
                pattern = Chem.MolFromSmarts(smarts_info["pattern"])
                if pattern:
                    found = mol.GetSubstructMatches(pattern)
                    if found:
                        matches.append({
                            "pattern": smarts_info["pattern"],
                            "description": smarts_info["desc"],
                            "count": len(found),
                        })

            alerts[tox_type] = {
                "name": info["name"],
                "description": info["description"],
                "risk": info["risk"],
                "alerts_found": len(matches),
                "matches": matches,
                "risk_level": "HIGH" if len(matches) >= 2 else "MEDIUM" if len(matches) == 1 else "LOW",
            }

        # Общая оценка
        total_alerts = sum(a["alerts_found"] for a in alerts.values())
        overall_risk = "HIGH" if total_alerts >= 4 else "MEDIUM" if total_alerts >= 2 else "LOW"

        return {
            "smiles": smiles,
            "overall_risk": overall_risk,
            "total_alerts": total_alerts,
            "assessments": alerts,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Ретросинтез
# ─────────────────────────────────────────────────────────────────────────────


# Типичные ретросинтетические разрывы
RETROSYNTHETIC_DISCONNECTIONS = {
    "ester_hydrolysis": {
        "name": "Гидролиз сложного эфира",
        "smarts": "[CX3](=O)[OX2H0][#6]>>[CX3](=O)[OX2H1].[OX2H0][#6]",
        "description": "Разрыв сложноэфирной связи на карбоновую кислоту и спирт",
    },
    "amide_hydrolysis": {
        "name": "Гидролиз амида",
        "smarts": "[CX3](=O)[NX3]>>[CX3](=O)[OX2H1].[NX3]",
        "description": "Разрыв амидной связи на карбоновую кислоту и амин",
    },
    "suzuki_coupling": {
        "name": "Сочетание Сузуки",
        "smarts": "[c,c][c,c]>>[c,c][B](O)O.[c,c][Cl,Br,I]",
        "description": "Разрыв биарильной связи на бороновую кислоту и арилгалогенид",
    },
    "reductive_amination": {
        "name": "Восстановительное аминирование",
        "smarts": "[CX3](=O)[NX3]>>[CX3](=O)[OX2H1].[NX3]",
        "description": "Разрыв связи C-N на карбонильное соединение и амин",
    },
    "michael_addition": {
        "name": "Присоединение Михаэля",
        "smarts": "[CX3](=O)[CX3]=[CX3]>>[CX3](=O)[CX3]=[CX3]",
        "description": "Разрыв α,β-ненасыщенного карбонильного соединения",
    },
    "wittig_reaction": {
        "name": "Реакция Виттига",
        "smarts": "[CX3]=[CX3]>>[CX3](=O).[CX3]=[PX3]",
        "description": "Разрыв двойной связи на карбонильное соединение и фосфониевый илид",
    },
    "grignard_addition": {
        "name": "Присоединение реактива Гриньяра",
        "smarts": "[CX3](=O)[#6]>>[CX3](=O)[OX2H1].[#6][Mg][Cl,Br]",
        "description": "Разрыв связи C-C рядом с карбонильной группой",
    },
    "friedel_crafts": {
        "name": "Реакция Фриделя-Крафтса",
        "smarts": "[c][CX3](=O)[Cl]>>[c].[CX3](=O)[Cl]",
        "description": "Разрыв связи арил-ацил на ароматическое соединение и ацилхлорид",
    },
}


@router.get("/retrosynthesis", summary="Ретросинтетический анализ")
def retrosynthetic_analysis(smiles: str = Query(..., description="SMILES целевой молекулы")):
    """
    Ретросинтетический анализ молекулы.
    Предлагает возможные разрывы связей для синтеза из коммерчески доступных прекурсоров.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = m.m
        suggestions = []

        for key, info in RETROSYNTHETIC_DISCONNECTIONS.items():
            rxn = AllChem.ReactionFromSmarts(info["smarts"])
            if rxn is None:
                continue

            # Прямая реакция (для проверки применимости)
            products = rxn.RunReactants((mol,))
            if products:
                # Проверяем, что продукты осмысленны
                valid_products = []
                for product_set in products:
                    for product in product_set:
                        try:
                            Chem.SanitizeMol(product)
                            smi = Chem.MolToSmiles(product, canonical=True)
                            if smi and "." not in smi:  # Только одиночные молекулы
                                valid_products.append(smi)
                        except Exception:
                            continue

                if valid_products:
                    suggestions.append({
                        "key": key,
                        "name": info["name"],
                        "description": info["description"],
                        "smarts": info["smarts"],
                        "possible_precursors": list(set(valid_products))[:5],  # До 5 вариантов
                    })

        return {
            "smiles": smiles,
            "suggestions_count": len(suggestions),
            "suggestions": suggestions,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Конформеры
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/conformers", summary="Генерация конформеров молекулы")
def generate_conformers(
    smiles: str = Query(..., description="SMILES молекулы"),
    num_conformers: int = Query(10, ge=1, le=50, description="Количество конформеров"),
):
    """
    Генерирует несколько энергетически выгодных конформеров молекулы.
    Использует ETKDG для генерации и MMFF для оптимизации.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = Chem.AddHs(m.m)

        # Генерация конформеров с ETKDG
        params = AllChem.ETKDGv3()
        params.numThreads = 0  # Использовать все доступные потоки
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

        if len(cids) == 0:
            raise ValueError("Failed to generate conformers")

        # Оптимизация каждого конформера
        results = []
        for cid in cids:
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=cid)
                energy = AllChem.MMFFGetMoleculeForceField(mol, confId=cid).CalcEnergy()
                results.append({
                    "conf_id": int(cid),
                    "energy": round(energy, 4),
                })
            except Exception:
                continue

        # Сортировка по энергии
        results.sort(key=lambda x: x["energy"])

        # Генерация SDF для каждого конформера
        conformers_sdf = []
        for i, r in enumerate(results):
            conf_mol = Chem.Mol(mol)
            conf_mol.RemoveAllConformers()
            conf_mol.AddConformer(mol.GetConformer(r["conf_id"]))
            sdf_block = Chem.MolToMolBlock(conf_mol)
            conformers_sdf.append({
                "rank": i + 1,
                "conf_id": r["conf_id"],
                "energy": r["energy"],
                "sdf": sdf_block,
            })

        return {
            "smiles": smiles,
            "total_conformers": len(results),
            "conformers": conformers_sdf,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Экспорт в SVG
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/export_svg", summary="Экспорт 2D-структуры в SVG")
def export_svg(
    smiles: str = Query(..., description="SMILES молекулы"),
    width: int = Query(400, ge=100, le=2000, description="Ширина в пикселях"),
    height: int = Query(400, ge=100, le=2000, description="Высота в пикселях"),
):
    """
    Экспортирует 2D-структуру молекулы в формате SVG.
    Подходит для вставки в научные публикации.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = m.m

        # Генерация 2D-координат
        AllChem.Compute2DCoords(mol)

        # Рендеринг в SVG
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        return Response(content=svg, media_type="image/svg+xml")
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Масс-спектрометрия (изотопное распределение)
# ─────────────────────────────────────────────────────────────────────────────
@router.get("/ms", summary="Масс-спектрометрия (изотопное распределение)")
def mass_spectrometry(smiles: str = Query(..., description="SMILES молекулы")):
    """
    Расчёт ожидаемого масс-спектра и изотопного распределения через pyOpenMS.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol = m.m
        exact_mass = Descriptors.ExactMolWt(mol)
        formula_str = rdMolDescriptors.CalcMolFormula(mol)

        import pyopenms as oms

        # Создаем формулу напрямую из корректной строки RDKit
        formula = oms.EmpiricalFormula(formula_str)

        # Генератор изотопных паттернов (4 пика)
        generator = oms.CoarseIsotopePatternGenerator(4)
        isotope_distribution = formula.getIsotopeDistribution(generator)

        # Извлекаем пики (переводим интенсивность в проценты для наглядности)
        peaks = []
        for peak in isotope_distribution.getContainer():
            peaks.append({
                "m_z": round(peak.getMZ(), 4),
                "intensity_percent": round(peak.getIntensity() * 100, 4)
            })

        return {
            "smiles": smiles,
            "exact_mass": round(exact_mass, 6),
            "molecular_formula": formula_str,
            "isotope_pattern": peaks,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

def _calculate_isotope_pattern(atom_counts, exact_mass):
    """Расчёт изотопного распределения с использованием pyOpenMS."""
    try:
        import pyopenms as oms

        # Формируем формулу строкой
        formula_str = ""
        for elem, count in sorted(atom_counts.items()):
            if count > 0:
                formula_str += elem
                if count > 1:
                    formula_str += str(count)

        # Создаём формулу
        formula = oms.EmpiricalFormula(formula_str)

        # Генератор изотопных паттернов (4 пика)
        generator = oms.CoarseIsotopePatternGenerator(4)
        isotope_distribution = formula.getIsotopeDistribution(generator)

        # Извлекаем пики
        peaks = []
        for peak in isotope_distribution.getContainer():
            peaks.append((round(peak.getMZ(), 4), round(peak.getIntensity(), 6)))

        # Сортировка по интенсивности
        peaks.sort(key=lambda x: x[1], reverse=True)

        # Нормализация к максимуму
        if peaks:
            max_intensity = peaks[0][1]
            peaks = [(m, round(i / max_intensity, 4)) for m, i in peaks]

        return peaks

    except ImportError:
        # Fallback на упрощённый расчёт, если pyOpenMS недоступен
        return _calculate_isotope_pattern_fallback(atom_counts, exact_mass)


def _calculate_isotope_pattern_fallback(atom_counts, exact_mass):
    """Упрощённый расчёт изотопного распределения (fallback)."""
    isotopes = {
        "C": [(12.0, 0.989), (13.003, 0.011)],
        "H": [(1.008, 0.9998), (2.014, 0.0002)],
        "N": [(14.003, 0.996), (15.000, 0.004)],
        "O": [(15.995, 0.998), (17.999, 0.002)],
        "S": [(31.972, 0.95), (33.968, 0.04)],
        "Cl": [(34.969, 0.75), (36.966, 0.25)],
        "Br": [(78.918, 0.50), (80.916, 0.50)],
    }

    peaks = [(round(exact_mass, 4), 1.0)]

    for elem, count in atom_counts.items():
        if elem not in isotopes:
            continue
        elem_isotopes = isotopes[elem]
        if len(elem_isotopes) < 2:
            continue

        new_peaks = []
        for mass, intensity in peaks:
            new_peaks.append((mass, intensity * elem_isotopes[0][1] ** count))
            if count >= 1:
                delta_mass = elem_isotopes[1][0] - elem_isotopes[0][0]
                new_peaks.append((
                    round(mass + delta_mass * count, 4),
                    intensity * elem_isotopes[1][1] * count
                ))
        peaks = new_peaks

    peaks.sort(key=lambda x: x[1], reverse=True)

    max_intensity = peaks[0][1] if peaks else 1.0
    normalized = [(m, round(i / max_intensity, 4)) for m, i in peaks[:10]]

    return normalized
