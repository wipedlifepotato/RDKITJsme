"""
Enhanced RDKit API Router — расширенный функционал хемоинформатики
"""

from fastapi import APIRouter, HTTPException, Query, Depends, Response
from fastapi.responses import HTMLResponse
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from pydantic import BaseModel, Field
from typing import List, Optional
import json
import math
import html as html_lib
import urllib.parse
import urllib.request

from database import get_db, CachedName
from Molecule import Molecule
from i18n import get_text, get_lang

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Recap
from rdkit.Chem import QED, Draw, inchi
from rdkit.Chem import Lipinski, Crippen, MolSurf, GraphDescriptors
from rdkit.Chem import rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolTransforms

router = APIRouter(prefix="/rdkit/api", tags=["RDKit Enhanced API"])

PROXY_ADDRESS = "127.0.0.1:9050"


@router.get("/", include_in_schema=False)
async def rdkit_index(lang: str = Depends(get_lang)):
    """Индексная страница RDKit API со списком всех эндпоинтов."""
    endpoints = {
        get_text("rdkit_index_adme", lang): {
            "GET /rdkit/api/adme": get_text("ep_adme_desc", lang),
        },
        get_text("rdkit_index_analysis", lang): {
            "GET /rdkit/api/functional_groups": get_text("ep_functional_groups_desc", lang),
            "GET /rdkit/api/validate": get_text("ep_validate_desc", lang),
            "GET /rdkit/api/pka": get_text("ep_pka_desc", lang),
            "GET /rdkit/api/descriptors": get_text("ep_descriptors_desc", lang),
        },
        get_text("rdkit_index_reactions", lang): {
            "GET /rdkit/api/reaction": get_text("ep_reaction_desc", lang),
            "GET /rdkit/api/by_inchikey": get_text("ep_by_inchikey_desc", lang),
        },
        get_text("rdkit_index_3d", lang): {
            "GET /rdkit/api/3d": get_text("ep_3d_desc", lang),
            "GET /rdkit/api/render_highlighted": get_text("ep_render_highlighted_desc", lang),
        },
        get_text("rdkit_index_batch", lang): {
            "POST /rdkit/api/batch": get_text("ep_batch_desc", lang),
            "GET /rdkit/api/library": get_text("ep_library_desc", lang),
            "POST /rdkit/api/library/save": get_text("ep_library_save_desc", lang),
            "DELETE /rdkit/api/library/{smiles}": get_text("ep_library_delete_desc", lang),
        },
        get_text("rdkit_index_comparison", lang): {
            "POST /rdkit/api/compare": get_text("ep_compare_desc", lang),
        },
    }
    return {"endpoints": endpoints, "total": sum(len(v) for v in endpoints.values())}


def _optimize_and_energy(mol, conf_id):
    """Оптимизация конформера: MMFF, а если параметров нет — UFF. Возвращает энергию (ккал/моль)."""
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=2000)
        props = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    else:
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=2000)
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
    return ff.CalcEnergy()


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
def get_adme_properties(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
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
            raise ValueError(get_text("invalid_smiles", lang))

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
            qed_interp = get_text("qed_excellent", lang)
        elif qed_score >= 0.5:
            qed_interp = get_text("qed_good", lang)
        elif qed_score >= 0.3:
            qed_interp = get_text("qed_fair", lang)
        else:
            qed_interp = get_text("qed_poor", lang)

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
    "nitro": {"smarts": "[$([NX3+](=O)[O-]),$([NX3](=O)=O)]", "name": "Nitro (-NO2)"},
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
def find_functional_groups(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Определяет наличие функциональных групп в молекуле по SMARTS-паттернам.
    Возвращает список найденных групп с количеством совпадений.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

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
                        "name": get_text(f"fg_{key}", lang),
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
def validate_smiles(
    smiles: str = Query(..., description=get_text("param_smiles_to_validate")),
    lang: str = Depends(get_lang),
):
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
        result["error_message"] = get_text("val_empty", lang)
        return result

    smiles = smiles.strip()

    # Проверка на недопустимые символы
    allowed_chars = set(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]=#@+-./\\:%*"
    )
    invalid_chars = set()
    for i, ch in enumerate(smiles):
        if ch not in allowed_chars:
            invalid_chars.add((i, ch))
    if invalid_chars:
        result["error_type"] = "invalid_characters"
        result["error_message"] = f"{get_text('val_invalid_characters', lang)}: {invalid_chars}"
        result["error_position"] = min(invalid_chars)[0]
        return result

    # Проверка баланса скобок
    open_paren = smiles.count("(")
    close_paren = smiles.count(")")
    if open_paren != close_paren:
        result["error_type"] = "unbalanced_parentheses"
        result["error_message"] = (
            f"{get_text('val_unbalanced_parentheses', lang)}: {open_paren} open vs {close_paren} close"
        )
        return result

    open_bracket = smiles.count("[")
    close_bracket = smiles.count("]")
    if open_bracket != close_bracket:
        result["error_type"] = "unbalanced_brackets"
        result["error_message"] = (
            f"{get_text('val_unbalanced_brackets', lang)}: {open_bracket} open vs {close_bracket} close"
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
                result["error_message"] = get_text("val_syntax_error", lang)
        except Exception:
            result["error_type"] = "syntax_error"
            result["error_message"] = get_text("val_syntax_error", lang)
        return result

    # Валидный SMILES
    result["valid"] = True
    result["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)

    # Предупреждения
    if mol.GetNumAtoms() > 100:
        result["warnings"].append(get_text("val_large_molecule", lang))
    if Descriptors.MolWt(mol) > 1000:
        result["warnings"].append(get_text("val_high_mw", lang))
    if Descriptors.MolLogP(mol) > 7:
        result["warnings"].append(get_text("val_high_logp", lang))

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Реакции (SMIRKS)
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/reaction", summary="Применение реакции (SMIRKS)")
def apply_reaction(
    smiles: str = Query(..., description=get_text("param_smiles_reagent")),
    reaction_smarts: str = Query(..., description=get_text("param_reaction_smarts")),
    lang: str = Depends(get_lang),
):
    """
    Применяет реакцию (SMIRKS) к молекуле.
    Возвращает список продуктов реакции.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(get_text("invalid_smiles", lang))

        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            raise ValueError(get_text("invalid_smirks_reaction", lang))

        products = rxn.RunReactants((mol,))

        if not products:
            return {
                "smiles": smiles,
                "reaction_smarts": reaction_smarts,
                "products": [],
                "product_count": 0,
                "message": get_text("no_reaction_products", lang),
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
def estimate_pka(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Оценивает кислотные/основные центры молекулы.
    Использует эвристический подход на основе функциональных групп.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

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
                        "type_name": get_text("pka_carboxylic_acid", lang),
                        "estimated_pka": 4.5,
                        "strength": get_text("pka_acidic", lang),
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
                        "type_name": get_text("pka_phenol", lang),
                        "estimated_pka": 10.0,
                        "strength": get_text("pka_acidic", lang),
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
                        "type_name": get_text("pka_aliphatic_amine", lang),
                        "estimated_pka": 10.0,
                        "strength": get_text("pka_basic", lang),
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
                        "type_name": get_text("pka_aniline", lang),
                        "estimated_pka": 4.6,
                        "strength": get_text("pka_basic", lang),
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
                        "type_name": get_text("pka_thiol", lang),
                        "estimated_pka": 9.0,
                        "strength": get_text("pka_acidic", lang),
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
                        "type_name": get_text("pka_sulfonic_acid", lang),
                        "estimated_pka": -2.0,
                        "strength": get_text("pka_strong_acid", lang),
                    }
                )

        return {
            "smiles": smiles,
            "pka_sites": pka_sites,
            "total_sites": len(pka_sites),
            "acidic_sites": len(
                [s for s in pka_sites if s["strength"] in (get_text("pka_acidic", lang), get_text("pka_strong_acid", lang))]
            ),
            "basic_sites": len([s for s in pka_sites if s["strength"] == get_text("pka_basic", lang)]),
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
def batch_analyze(request: BatchRequest, lang: str = Depends(get_lang)):
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
                    {"smiles": smiles, "valid": False, "error": get_text("invalid_smiles", lang)}
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
def get_by_inchikey(
    inchikey: str = Query(..., description=get_text("param_inchikey")),
    lang: str = Depends(get_lang),
):
    """
    Получает свойства молекулы по InChIKey.
    Использует Cactus/NCI resolver для конвертации.
    """
    try:
        # Проверка формата InChIKey
        if len(inchikey) != 27 or inchikey.count("-") != 2:
            raise ValueError(get_text("invalid_inchikey_format", lang))

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
            raise ValueError(get_text("could_not_resolve_inchikey", lang))

        # Получаем свойства
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("resolved_smiles_invalid", lang))

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
def get_all_descriptors(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Возвращает все доступные дескрипторы RDKit для молекулы.
    Включает 200+ молекулярных дескрипторов.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol = m.m

        # Вычисляем все дескрипторы
        descriptors = {}
        for name, func in Descriptors.descList:
            try:
                value = func(mol)
                if isinstance(value, (int, float)):
                    if isinstance(value, float) and not math.isfinite(value):
                        descriptors[name] = None
                    else:
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

        extra = {
            k: (None if isinstance(v, float) and not math.isfinite(v) else v)
            for k, v in extra.items()
        }

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
    smiles: str = Query(..., description=get_text("param_smiles")),
    format: str = Query("sdf", description=get_text("param_format")),
    lang: str = Depends(get_lang),
):
    """
    Генерирует 3D-структуру молекулы.
    Возвращает в формате SDF, MOL или PDB.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol = Chem.AddHs(m.m)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) == -1:
            raise ValueError(get_text("failed_to_embed_molecule", lang))
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)

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
            raise ValueError(f"{get_text('unsupported_format', lang)}: {format}")

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
def save_to_library(
    entry: MoleculeEntry,
    lang: str = Depends(get_lang),
    db: Session = Depends(get_db),
):
    """Сохраняет молекулу в локальную библиотеку (SQLite)."""
    try:
        m = Molecule(entry.smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        canonical = Chem.MolToSmiles(m.m, canonical=True)

        # Проверяем, есть ли уже
        existing = db.query(CachedName).filter(CachedName.smiles == canonical).first()
        if existing:
            if entry.name:
                existing.name = entry.name
            db.commit()
            return {"status": get_text("lib_updated", lang), "smiles": canonical}

        new_entry = CachedName(smiles=canonical, name=entry.name or "")
        db.add(new_entry)
        db.commit()
        return {"status": get_text("lib_saved", lang), "smiles": canonical}
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
def delete_from_library(
    smiles: str,
    lang: str = Depends(get_lang),
    db: Session = Depends(get_db),
):
    """Удаляет молекулу из библиотеки по SMILES."""
    try:
        entry = db.query(CachedName).filter(CachedName.smiles == smiles).first()
        if not entry:
            raise HTTPException(status_code=404, detail=get_text("molecule_not_found_in_library", lang))
        db.delete(entry)
        db.commit()
        return {"status": get_text("lib_deleted", lang), "smiles": smiles}
    except HTTPException:
        raise
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Сравнение нескольких молекул
# ─────────────────────────────────────────────────────────────────────────────


@router.post("/compare", summary="Сравнение нескольких молекул")
def compare_molecules(
    smiles_list: List[str] = Query(..., description=get_text("param_smiles_list")),
    lang: str = Depends(get_lang),
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
            raise ValueError(get_text("need_at_least_two_smiles", lang))

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
    smiles: str = Query(..., description=get_text("param_smiles")),
    highlight_groups: Optional[List[str]] = Query(
        None, description=get_text("param_highlight_groups")
    ),
    lang: str = Depends(get_lang),
):
    """
    Рендерит молекулу с подсветкой функциональных групп.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

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
def estimate_nmr(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Эвристическая оценка химических сдвигов ЯМР с учетом химического окружения.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

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
                            "description": get_text("nmr_cooh", lang),
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
                    desc = get_text("nmr_ch3_near_carbonyl", lang)
                else:
                    shift_range = [0.8, 1.2]
                    desc = get_text("nmr_ch3_aliphatic", lang)

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
                c_shift = _estimate_13c_shift_improved(mol, atom, lang)
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


def _estimate_13c_shift_improved(mol, c_atom, lang="ru"):
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
                return {"range": [165, 185], "desc": get_text("nmr_carbonyl_acid_ester", lang)}
            return {"range": [190, 220], "desc": get_text("nmr_carbonyl_ketone_aldehyde", lang)}

    # Ароматические углероды
    if c_atom.GetIsAromatic():
        return {"range": [110, 160], "desc": get_text("nmr_aromatic_c", lang)}

    # Алифатические углероды по количеству водородов
    num_h = sum(1 for n in c_atom.GetNeighbors() if n.GetAtomicNum() == 1)
    if num_h >= 3:
        return {"range": [10, 30], "desc": get_text("nmr_methyl", lang)}
    elif num_h == 2:
        return {"range": [15, 45], "desc": get_text("nmr_methylene", lang)}
    else:
        return {"range": [25, 60], "desc": get_text("nmr_ch_methine_quaternary", lang)}

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


@router.get("/toxicity", summary="Оценка токсичности через каталоги RDKit (Brenk / PAINS)")
def assess_toxicity(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Оценка потенциальной токсичности молекулы через FilterCatalog RDKit:
    - Brenk filters (нежелательные фрагменты для лекарств)
    - PAINS (пан-ассайные интерферирующие вещества)
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))
        mol = m.m

        # Подключаем каталоги фильтров (Brenk и PAINS)
        from rdkit.Chem import FilterCatalog
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)

        matches = catalog.GetMatches(mol)
        alerts_found = []
        for entry in matches:
            alerts_found.append({
                "name": entry.GetDescription(),
                "heading": entry.GetHeading() if hasattr(entry, 'GetHeading') else get_text("tox_alert", lang)
            })

        if len(alerts_found) >= 2:
            overall_risk = get_text("tox_risk_high", lang)
        elif len(alerts_found) == 1:
            overall_risk = get_text("tox_risk_medium", lang)
        else:
            overall_risk = get_text("tox_risk_low", lang)

        return {
            "smiles": smiles,
            "overall_risk": overall_risk,
            "total_alerts": len(alerts_found),
            "alerts": alerts_found
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Ретросинтез (RDKit RECAP)
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/retrosynthesis", summary="Ретросинтетический анализ (RECAP)")
def retrosynthetic_analysis(
    smiles: str = Query(..., description=get_text("param_smiles_target")),
    lang: str = Depends(get_lang),
):
    """
    Ретросинтетический анализ молекулы алгоритмом RECAP
    (Retrosynthetic Combinatorial Analysis Procedure).

    Разрывает связи по стандартным правилам (амиды, сложные эфиры,
    простые эфиры, мочевины, сульфонамиды и т.д.) и возвращает конечные
    фрагменты-прекурсоры — только листья иерархического дерева декомпозиции.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol = m.m

        # RecapDecompose возвращает корень иерархического дерева фрагментов
        # (или None, если молекулу не удалось разобрать).
        recap_tree = Recap.RecapDecompose(mol)
        if recap_tree is None:
            return {
                "smiles": smiles,
                "decomposable": False,
                "fragments_count": 0,
                "fragments": [],
            }

        # GetLeaves() — только конечные фрагменты (листья дерева).
        # GetAllChildren() вернул бы ВСЕ узлы, включая промежуточные.
        leaves = recap_tree.GetLeaves()

        fragments = []
        seen = set()
        for leaf_smiles in leaves:
            frag_mol = Chem.MolFromSmiles(leaf_smiles)
            if frag_mol is None:
                continue
            # Канонизируем SMILES для единообразия и дедупликации
            canon_smi = Chem.MolToSmiles(frag_mol)
            if canon_smi in seen:
                continue
            seen.add(canon_smi)
            fragments.append(
                {
                    "smiles": canon_smi,
                    "formula": rdMolDescriptors.CalcMolFormula(frag_mol),
                    "molecular_weight": round(Descriptors.MolWt(frag_mol), 4),
                    "attachment_points": canon_smi.count("*"),
                }
            )

        return {
            "smiles": smiles,
            "decomposable": bool(fragments),
            "fragments_count": len(fragments),
            "fragments": fragments,
            "note": get_text("retro_note", lang),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Конформеры
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/conformers", summary="Генерация конформеров молекулы")
def generate_conformers(
    smiles: str = Query(..., description=get_text("param_smiles")),
    num_conformers: int = Query(10, ge=1, le=50, description=get_text("param_num_conformers")),
    lang: str = Depends(get_lang),
):
    """
    Генерирует несколько энергетически выгодных конформеров молекулы.
    Использует ETKDG для генерации и MMFF для оптимизации.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol = Chem.AddHs(m.m)

        # Генерация конформеров с ETKDG
        params = AllChem.ETKDGv3()
        params.numThreads = 0  # Использовать все доступные потоки
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

        if len(cids) == 0:
            raise ValueError(get_text("failed_to_generate_conformers", lang))

        # Оптимизация каждого конформера
        results = []
        for cid in cids:
            try:
                energy = _optimize_and_energy(mol, cid)
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
    smiles: str = Query(..., description=get_text("param_smiles")),
    width: int = Query(400, ge=100, le=2000, description=get_text("param_width")),
    height: int = Query(400, ge=100, le=2000, description=get_text("param_height")),
    lang: str = Depends(get_lang),
):
    """
    Экспортирует 2D-структуру молекулы в формате SVG.
    Подходит для вставки в научные публикации.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

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
def mass_spectrometry(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Расчёт ожидаемого масс-спектра и изотопного распределения через pyOpenMS.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

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


# ─────────────────────────────────────────────────────────────────────────────
# Тавтомеры
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/tautomers", summary="Генерация таутомеров")
def enumerate_tautomers(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Перечисляет все возможные таутомеры молекулы.

    Использует алгоритм TautomerEnumerator из RDKit для генерации
    всех уникальных таутомерных форм молекулы.

    **Параметры:**
    - `smiles` — строка SMILES молекулы
    - `lang` — язык ответа (ru, en, uk, es)

    **Возвращает:**
    - `canonical_tautomer` — канонический (наиболее стабильный) таутомер
    - `tautomers` — список всех уникальных таутомеров в формате SMILES
    - `count` — количество найденных таутомеров

    **Пример запроса:**
    ```
    GET /rdkit/api/tautomers?smiles=CC(=O)C&lang=ru
    ```

    **Пример ответа:**
    ```json
    {
        "smiles": "CC(=O)C",
        "canonical_tautomer": "CC(=O)C",
        "tautomers": ["CC(=O)C", "C=C(O)C"],
        "count": 2
    }
    ```
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        enumerator = rdMolStandardize.TautomerEnumerator()
        canonical = enumerator.Canonicalize(m.m)
        tautomers = enumerator.Enumerate(m.m)

        canonical_smiles = Chem.MolToSmiles(canonical, canonical=True)
        all_tautomers = [Chem.MolToSmiles(t, canonical=True) for t in tautomers]

        return {
            "smiles": smiles,
            "canonical_tautomer": canonical_smiles,
            "tautomers": all_tautomers,
            "count": len(all_tautomers),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Bemis-Murcko Scaffold
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/scaffold", summary="Выделение скаффолда Bemis-Murcko")
def get_scaffold(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Выделяет каркас Bemis-Murcko из молекулы.

    Алгоритм Bemis-Murcko удаляет боковые цепи, оставляя только
    циклический скелет с связями между циклами.

    **Параметры:**
    - `smiles` — строка SMILES молекулы
    - `lang` — язык ответа (ru, en, uk, es)

    **Возвращает:**
    - `scaffold_smiles` — SMILES скаффолда
    - `scaffold_svg` — SVG-изображение скаффолда

    **Пример запроса:**
    ```
    GET /rdkit/api/scaffold?smiles=CC(=O)Oc1ccccc1C(=O)O&lang=ru
    ```

    **Пример ответа:**
    ```json
    {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "scaffold_smiles": "O=C(O)c1ccccc1",
        "scaffold_svg": "<svg>...</svg>"
    }
    ```
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        scaffold = MurckoScaffold.GetScaffoldForMol(m.m)
        scaffold_smiles = Chem.MolToSmiles(scaffold, canonical=True)

        # Генерация SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(scaffold)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        return {
            "smiles": smiles,
            "scaffold_smiles": scaffold_smiles,
            "scaffold_svg": svg,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Butina Clustering
# ─────────────────────────────────────────────────────────────────────────────


class ClusterRequest(BaseModel):
    """Модель запроса для кластеризации молекул."""

    smiles_list: List[str] = Field(..., description="Список SMILES строк (минимум 2)")
    cutoff: float = Field(0.7, ge=0.0, le=1.0, description="Порог сходства Tanimoto (0.0–1.0)")


@router.post("/cluster", summary="Кластеризация молекул (Butina)")
def cluster_molecules(
    request: ClusterRequest,
    lang: str = Depends(get_lang),
):
    """
    Кластеризует молекулы алгоритмом Butina.

    Использует Morgan Fingerprints (радиус 2, 2048 бит) и матрицу
    сходства Tanimoto для группировки молекул по структурному сходству.

    **Параметры запроса (JSON):**
    - `smiles_list` — список SMILES строк (минимум 2)
    - `cutoff` — порог сходства Tanimoto (0.0–1.0, по умолчанию 0.7)

    **Возвращает:**
    - `total_molecules` — количество валидных молекул
    - `cutoff` — использованный порог
    - `cluster_count` — количество кластеров
    - `clusters` — список кластеров с SMILES

    **Пример запроса:**
    ```json
    POST /rdkit/api/cluster
    {
        "smiles_list": ["CCO", "CCCO", "c1ccccc1", "c1ccccc1O"],
        "cutoff": 0.7
    }
    ```

    **Пример ответа:**
    ```json
    {
        "total_molecules": 4,
        "cutoff": 0.7,
        "cluster_count": 2,
        "clusters": [
            {"cluster_id": 0, "size": 2, "smiles": ["CCO", "CCCO"]},
            {"cluster_id": 1, "size": 2, "smiles": ["c1ccccc1", "c1ccccc1O"]}
        ]
    }
    ```
    """
    try:
        mols = []
        valid_smiles = []
        for smiles in request.smiles_list:
            m = Molecule(smiles)
            if m.m is not None:
                mols.append(m.m)
                valid_smiles.append(smiles)

        if len(mols) < 2:
            raise ValueError(get_text("need_at_least_two_smiles", lang))

        fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]

        # Вычисление матрицы расстояний (1 - Tanimoto)
        dists = []
        for i in range(1, len(fps)):
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
            dists.extend([1 - x for x in sims])

        clusters = Butina.ClusterData(dists, len(fps), 1.0 - request.cutoff, isDistData=True)  # cutoff — сходство, Butina ждёт расстояние

        result = []
        for cluster in clusters:
            result.append(
                {
                    "cluster_id": len(result),
                    "size": len(cluster),
                    "smiles": [valid_smiles[i] for i in cluster],
                }
            )

        return {
            "total_molecules": len(valid_smiles),
            "cutoff": request.cutoff,
            "cluster_count": len(clusters),
            "clusters": result,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Стереоизомеры
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/stereoisomers", summary="Генерация стереоизомеров")
def enumerate_stereoisomers(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Перечисляет все возможные стереоизомеры молекулы.

    Использует алгоритм EnumerateStereoisomers из RDKit для генерации
    всех уникальных стереоизомеров с неопределённой стереохимией.

    **Параметры:**
    - `smiles` — строка SMILES молекулы
    - `lang` — язык ответа (ru, en, uk, es)

    **Возвращает:**
    - `stereoisomers` — список всех стереоизомеров в формате SMILES
    - `count` — количество найденных стереоизомеров
    - `chiral_centers_count` — количество хиральных центров (N)
    - `chiral_centers` — список хиральных центров
    - `max_theoretical_stereoisomers` — максимальное теоретическое количество (2^N)
    - `formula` — строка с формулой (например, "2^2 = 4")

    **Проверка по формуле:**
    Количество стереоизомеров = 2^N, где N — количество хиральных центров.
    Например, для молекулы с 2 хиральными центрами: 2^2 = 4 стереоизомера.

    **Пример запроса:**
    ```
    GET /rdkit/api/stereoisomers?smiles=CC(O)C&lang=ru
    ```

    **Пример ответа:**
    ```json
    {
        "smiles": "CC(O)C",
        "stereoisomers": ["C[C@H](O)C", "C[C@@H](O)C"],
        "count": 2,
        "chiral_centers_count": 1,
        "chiral_centers": [[1, "?")],
        "max_theoretical_stereoisomers": 2,
        "formula": "2^1 = 2"
    }
    ```
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        # Находим хиральные центры
        chiral_centers = Chem.FindMolChiralCenters(m.m, includeUnassigned=True)
        num_chiral = len(chiral_centers)
        max_theoretical = 2 ** num_chiral if num_chiral > 0 else 1

        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True)
        isomers = list(EnumerateStereoisomers(m.m, options=opts))

        isomer_smiles = [
            Chem.MolToSmiles(iso, canonical=True, isomericSmiles=True) for iso in isomers
        ]

        return {
            "smiles": smiles,
            "stereoisomers": isomer_smiles,
            "count": len(isomer_smiles),
            "chiral_centers_count": num_chiral,
            "chiral_centers": chiral_centers,
            "max_theoretical_stereoisomers": max_theoretical,
            "formula": f"2^{num_chiral} = {max_theoretical}",
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# 3D-конформеры (интерактивная визуализация)
# ─────────────────────────────────────────────────────────────────────────────


@router.get("/3d_conformers", summary="3D-конформеры с интерактивной визуализацией")
def get_3d_conformers(
    smiles: str = Query(..., description=get_text("param_smiles")),
    num_conformers: int = Query(5, ge=1, le=20, description=get_text("param_num_conformers")),
    lang: str = Depends(get_lang),
):
    """
    Генерирует 3D-конформеры молекулы с интерактивной визуализацией.

    Возвращает HTML-страницу с 3Dmol.js визуализацией:
    - Вращение молекулы мышью
    - Отображение энергии каждого конформера
    - Переключение между конформерами
    - Автоматическое вращение (spin)

    **Параметры:**
    - `smiles` — строка SMILES молекулы
    - `num_conformers` — количество конформеров (1–20, по умолчанию 5)
    - `lang` — язык ответа (ru, en, uk, es)

    **Возвращает:**
    - HTML-страницу с интерактивной 3D-визуализацией

    **Пример запроса:**
    ```
    GET /rdkit/api/3d_conformers?smiles=CCO&num_conformers=5&lang=ru
    ```

    **Примечание:**
    Для просмотра откройте ответ в браузере. Страница использует 3Dmol.js
    для интерактивной визуализации молекулы.
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol = Chem.AddHs(m.m)

        # Генерация конформеров
        params = AllChem.ETKDGv3()
        params.numThreads = 0
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

        if len(cids) == 0:
            raise ValueError(get_text("failed_to_generate_conformers", lang))

        # Оптимизация и сбор данных
        conformers = []
        for cid in cids:
            try:
                energy = _optimize_and_energy(mol, cid)

                # Координаты атомов
                conf = mol.GetConformer(cid)
                atoms = []
                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    atoms.append(
                        {
                            "index": atom.GetIdx(),
                            "element": atom.GetSymbol(),
                            "x": round(pos.x, 4),
                            "y": round(pos.y, 4),
                            "z": round(pos.z, 4),
                        }
                    )

                # Связи
                bonds = []
                for bond in mol.GetBonds():
                    bonds.append(
                        {
                            "begin": bond.GetBeginAtomIdx(),
                            "end": bond.GetEndAtomIdx(),
                            "type": str(bond.GetBondType()),
                        }
                    )

                conformers.append(
                    {
                        "conf_id": int(cid),
                        "energy": round(energy, 4),
                        "atoms": atoms,
                        "bonds": bonds,
                    }
                )
            except Exception:
                continue

        # Сортировка по энергии
        conformers.sort(key=lambda x: x["energy"])
        if not conformers:
            raise ValueError(get_text("failed_to_generate_conformers", lang))

        conformers_json = json.dumps(conformers).replace("</", "<\\/")
        smiles_html = html_lib.escape(smiles)

        # Генерация HTML с 3Dmol.js
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>3D Conformers — {smiles_html}</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        #viewer {{ width: 600px; height: 500px; position: relative; }}
        .controls {{ margin: 10px 0; }}
        .info {{ margin: 10px 0; padding: 10px; background: #f0f0f0; border-radius: 5px; }}
        button {{ padding: 8px 16px; margin: 5px; cursor: pointer; }}
        button.active {{ background: #4CAF50; color: white; }}
    </style>
</head>
<body>
    <h2>3D-конформеры: {smiles_html}</h2>
    <div class="info">
        <strong>Конформер:</strong> <span id="conf-id">1</span> |
        <strong>Энергия:</strong> <span id="conf-energy">{conformers[0]["energy"] if conformers else 0}</span> ккал/моль |
        <strong>Всего:</strong> {len(conformers)}
    </div>
    <div class="controls">
        <button onclick="prevConformer()">← Предыдущий</button>
        <button onclick="nextConformer()">Следующий →</button>
        <button onclick="toggleSpin()">Вращение</button>
        <button onclick="resetView()">Сброс</button>
    </div>
    <div id="viewer"></div>

    <script>
        const conformers = {conformers_json};
        let currentConf = 0;
        let spinning = false;
        let viewer = null;

        function initViewer() {{
            let element = document.getElementById('viewer');
            let config = {{ backgroundColor: 'white' }};
            viewer = $3Dmol.createViewer(element, config);
            loadConformer(0);
            viewer.zoomTo();
            viewer.render();
        }}

        function loadConformer(idx) {{
            currentConf = idx;
            let conf = conformers[idx];
            viewer.clear();

            // Добавляем атомы
            conf.atoms.forEach(atom => {{
                viewer.addAtom({{
                    x: atom.x,
                    y: atom.y,
                    z: atom.z,
                    elem: atom.element
                }});
            }});

            // Добавляем связи
            conf.bonds.forEach(bond => {{
                viewer.addBond({{ atom1: bond.begin, atom2: bond.end }});
            }});

            // Обновляем информацию
            document.getElementById('conf-id').textContent = idx + 1;
            document.getElementById('conf-energy').textContent = conf.energy;

            viewer.setStyle({{}}, {{stick: {{}}}});
            viewer.zoomTo();
            viewer.render();
        }}

        function nextConformer() {{
            let next = (currentConf + 1) % conformers.length;
            loadConformer(next);
        }}

        function prevConformer() {{
            let prev = (currentConf - 1 + conformers.length) % conformers.length;
            loadConformer(prev);
        }}

        function toggleSpin() {{
            spinning = !spinning;
            viewer.spin(spinning ? 'y' : false);
        }}

        function resetView() {{
            viewer.zoomTo();
            viewer.render();
        }}

        initViewer();
    </script>
</body>
</html>"""

        return HTMLResponse(content=html_content)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ─────────────────────────────────────────────────────────────────────────────
# Анализ конформаций колец (Ring Puckering)
# ─────────────────────────────────────────────────────────────────────────────


def _calculate_cremer_pople(coords):
    """
    Вычисляет параметры складчатости Кремера-Попла для цикла.
    
    Для 6-членных циклов: theta, phi, Q
    Для 5-членных циклов: theta, phi, q
    
    coords: список (x, y, z) координат атомов цикла
    """
    import math
    import numpy as np

    n = len(coords)
    if n < 5:
        return None

    # Центроид цикла
    centroid = np.mean(coords, axis=0)
    
    # Векторы от центроида к атомам
    vecs = coords - centroid
    
    # Для 6-членных циклов
    if n == 6:
        # Вычисляем z-координаты относительно средней плоскости
        # Используем SVD для нахождения нормали к плоскости
        _, _, vh = np.linalg.svd(vecs)
        normal = vh[-1]  # Последняя строка — нормаль к плоскости
        
        z = np.dot(vecs, normal)
        
        # Параметры Кремера-Попла для 6-членного цикла
        # Q = sqrt(sum(z_i^2))
        Q = np.sqrt(np.sum(z**2))
        
        # theta и phi
        if Q < 1e-6:
            return {"Q": 0.0, "theta": 0.0, "phi": 0.0, "conformation": "planar"}
        
        # Вычисляем theta и phi по формулам
        # Для 6-членного цикла:
        # z_i = Q * cos(theta) * cos(phi + 2*pi*i/6) для i=0..5
        # Используем преобразование Фурье
        
        # Упрощённый подход: используем первые два гармонических компонента
        angles = np.array([2 * np.pi * i / 6 for i in range(6)])
        
        # Разложение по косинусам и синусам
        a1 = np.sum(z * np.cos(angles))
        b1 = np.sum(z * np.sin(angles))
        a2 = np.sum(z * np.cos(2 * angles))
        b2 = np.sum(z * np.sin(2 * angles))
        
        # Для 6-членного цикла основная гармоника — m=2
        # theta = arccos(a2 / (Q * sqrt(3/2)))
        # phi = atan2(b2, a2)
        
        # Нормализация
        denom = Q * np.sqrt(1.5)
        if denom > 1e-6:
            cos_theta = np.clip(a2 / denom, -1, 1)
            theta = np.arccos(cos_theta) * 180 / np.pi
            phi = np.arctan2(b2, a2) * 180 / np.pi
        else:
            theta = 0.0
            phi = 0.0
        
        # Классификация конформации
        conformation = _classify_6ring(theta, phi)
        
        return {
            "Q": round(float(Q), 4),
            "theta": round(float(theta), 2),
            "phi": round(float(phi), 2),
            "conformation": conformation,
        }
    
    # Для 5-членных циклов
    elif n == 5:
        _, _, vh = np.linalg.svd(vecs)
        normal = vh[-1]
        z = np.dot(vecs, normal)
        
        q = np.sqrt(np.sum(z**2))
        
        if q < 1e-6:
            return {"q": 0.0, "theta": 0.0, "phi": 0.0, "conformation": "planar"}
        
        angles = np.array([2 * np.pi * i / 5 for i in range(5)])
        
        a1 = np.sum(z * np.cos(angles))
        b1 = np.sum(z * np.sin(angles))
        a2 = np.sum(z * np.cos(2 * angles))
        b2 = np.sum(z * np.sin(2 * angles))
        
        # Для 5-членного цикла основная гармоника — m=2
        denom = q * np.sqrt(2)
        if denom > 1e-6:
            cos_theta = np.clip(a2 / denom, -1, 1)
            theta = np.arccos(cos_theta) * 180 / np.pi
            phi = np.arctan2(b2, a2) * 180 / np.pi
        else:
            theta = 0.0
            phi = 0.0
        
        conformation = _classify_5ring(theta, phi)
        
        return {
            "q": round(float(q), 4),
            "theta": round(float(theta), 2),
            "phi": round(float(phi), 2),
            "conformation": conformation,
        }
    
    return None


def _classify_6ring(theta, phi):
    """Классификация конформации 6-членного цикла."""
    if theta < 30 or theta > 150:
        return "chair"
    elif 30 <= theta <= 60 or 120 <= theta <= 150:
        return "half-chair"
    elif 60 < theta < 120:
        # Различаем boat и twist-boat по phi
        phi_mod = phi % 60
        if phi_mod < 15 or phi_mod > 45:
            return "boat"
        else:
            return "twist-boat"
    return "unknown"


def _classify_5ring(theta, phi):
    """Классификация конформации 5-членного цикла."""
    if theta < 30:
        return "envelope"
    elif 30 <= theta <= 90:
        return "half-envelope"
    else:
        return "twist"
    return "unknown"


def _cremer_pople(coords):
    """
    Параметры складчатости Кремера–Попла (1975) по 3D-координатам атомов цикла
    (атомы должны идти по порядку обхода цикла). Поддерживаются 5- и 6-членные циклы.
    Q/q — амплитуда в Å, theta/phi — в градусах.
    """
    R = np.asarray(coords, dtype=float)
    n = len(R)
    if n not in (5, 6):
        return None
    R = R - R.mean(axis=0)
    ang = 2.0 * np.pi * np.arange(n) / n
    r1 = (R * np.sin(ang)[:, None]).sum(axis=0)
    r2 = (R * np.cos(ang)[:, None]).sum(axis=0)
    nv = np.cross(r1, r2)
    norm = np.linalg.norm(nv)
    if norm < 1e-8:
        return None
    z = R @ (nv / norm)  # отклонения от средней плоскости
    amp = float(np.sqrt(np.sum(z ** 2)))

    j = np.arange(n)
    a2 = np.sqrt(2.0 / n) * np.sum(z * np.cos(2 * ang))
    b2 = -np.sqrt(2.0 / n) * np.sum(z * np.sin(2 * ang))
    q2 = float(np.hypot(a2, b2))
    phi = float(np.degrees(np.arctan2(b2, a2)) % 360.0) if q2 > 0.02 else 0.0  # для кресла φ не определён

    if n == 6:
        q3 = float(np.sqrt(1.0 / n) * np.sum(((-1.0) ** j) * z))
        theta = float(np.degrees(np.arctan2(q2, q3)))
        if amp < 0.05:
            theta, phi = 0.0, 0.0  # у плоского цикла углы не определены — не показываем численный шум
            label = "planar"
        elif theta <= 25 or theta >= 155:
            label = "chair"
        elif 65 <= theta <= 115:
            pm = phi % 60.0
            label = "boat" if (pm <= 15 or pm >= 45) else "twist-boat"
        else:
            label = "half-chair/envelope"
        return {"Q": round(amp, 4), "theta": round(theta, 2), "phi": round(phi, 2), "conformation": label}

    # n == 5: единственная мода m=2, θ не определён
    pm = phi % 36.0
    if amp < 0.05:
        phi = 0.0
        label = "planar"
    else:
        label = "envelope" if (pm <= 9 or pm >= 27) else "twist"
    return {"q": round(amp, 4), "theta": None, "phi": round(phi, 2), "conformation": label}


def _calculate_cremer_pople_from_dihedrals(dihedrals, n):
    """
    Вычисляет параметры Кремера-Попла из диэдральных углов.
    Упрощённая версия для классификации конформаций.
    """
    import math

    if n == 6:
        # Для 6-членных циклов
        # Сумма квадратов диэдральных углов
        sum_sq = sum(d**2 for d in dihedrals)
        Q = math.sqrt(sum_sq / n)

        # Среднее значение
        mean = sum(dihedrals) / n

        # Определяем тип конформации по характерным признакам
        # Chair: чередующиеся +60/-60
        # Boat: два угла ~0, остальные ±60
        # Twist-boat: промежуточные значения

        # Считаем количество углов близких к 0 (характерно для boat)
        near_zero = sum(1 for d in dihedrals if abs(d) < 15)

        # Считаем количество углов близких к ±60 (характерно для chair)
        near_sixty = sum(1 for d in dihedrals if abs(abs(d) - 60) < 15)

        if near_sixty >= 4:
            conformation = "chair"
        elif near_zero >= 2:
            conformation = "boat"
        else:
            conformation = "twist-boat"

        # Вычисляем theta и phi на основе диэдральных углов
        # Для chair: theta ~ 0 или 180, phi ~ 0
        # Для boat: theta ~ 90, phi ~ 0 или 60
        # Для twist-boat: theta ~ 90, phi ~ 30 или 90

        if conformation == "chair":
            theta = 0.0 if mean > 0 else 180.0
            phi = 0.0
        elif conformation == "boat":
            theta = 90.0
            phi = 0.0 if near_zero >= 2 else 60.0
        else:  # twist-boat
            theta = 90.0
            phi = 30.0

        return {
            "Q": round(Q, 4),
            "theta": round(theta, 2),
            "phi": round(phi, 2),
            "conformation": conformation,
        }

    elif n == 5:
        # Для 5-членных циклов
        sum_sq = sum(d**2 for d in dihedrals)
        q = math.sqrt(sum_sq / n)

        # Определяем тип конформации
        # Envelope: один атом выбит из плоскости (один угол близок к 0)
        # Twist: два атома выбиты в разные стороны
        near_zero = sum(1 for d in dihedrals if abs(d) < 20)
        max_angle = max(abs(d) for d in dihedrals)

        if near_zero >= 1 and max_angle > 30:
            conformation = "envelope"
        elif max_angle > 30:
            conformation = "twist"
        else:
            conformation = "planar"

        return {
            "q": round(q, 4),
            "theta": 0.0,
            "phi": 0.0,
            "conformation": conformation,
        }

    return None


def _get_ring_conformation_for_conf(mol, ring_atoms, conf_id):
    """
    Определяет конформацию цикла для конкретного конформера.
    Использует rdMolTransforms.GetDihedralDeg для корректного вычисления углов.

    ring_atoms: список индексов атомов цикла
    conf_id: ID конформера
    """
    if len(ring_atoms) < 5:
        return None

    conf = mol.GetConformer(conf_id)
    n = len(ring_atoms)

    # Диэдральные углы через rdMolTransforms
    dihedrals = []
    for i in range(n):
        i1 = ring_atoms[i]
        i2 = ring_atoms[(i + 1) % n]
        i3 = ring_atoms[(i + 2) % n]
        i4 = ring_atoms[(i + 3) % n]
        angle = rdMolTransforms.GetDihedralDeg(conf, i1, i2, i3, i4)
        dihedrals.append(round(float(angle), 2))

    # Параметры Кремера-Попла
    coords = []
    for idx in ring_atoms:
        pos = conf.GetAtomPosition(idx)
        coords.append([pos.x, pos.y, pos.z])
    cp_params = _cremer_pople(coords)
    if cp_params is None:  # циклы других размеров — грубая эвристика по диэдрам
        cp_params = _calculate_cremer_pople_from_dihedrals(dihedrals, n)

    return {
        "conf_id": int(conf_id),
        "dihedral_angles": dihedrals,
        "cremer_pople": cp_params,
    }


@router.get("/ring_conformation", summary="Анализ конформаций колец")
def analyze_ring_conformation(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    """
    Анализирует конформации всех циклов в молекуле.

    Определяет тип конформации цикла (Chair, Boat, Twist-boat, Envelope и т.д.)
    и вычисляет параметры складчатости Кремера-Попла (θ, φ, Q).

    **Параметры:**
    - `smiles` — строка SMILES молекулы
    - `lang` — язык ответа (ru, en, uk, es)

    **Возвращает:**
    - `total_rings` — количество циклов
    - `rings` — список циклов с параметрами:
      - `ring_size` — размер цикла
      - `atom_indices` — индексы атомов цикла
      - `dihedral_angles` — диэдральные углы
      - `cremer_pople` — параметры Кремера-Попла:
        - `Q` — амплитуда складчатости (для 6-членных)
        - `q` — амплитуда складчатости (для 5-членных)
        - `theta` — полярный угол (°)
        - `phi` — азимутальный угол (°)
        - `conformation` — тип конформации

    **Пример запроса:**
    ```
    GET /rdkit/api/ring_conformation?smiles=C1CCCCC1&lang=ru
    ```

    **Пример ответа:**
    ```json
    {
        "smiles": "C1CCCCC1",
        "total_rings": 1,
        "rings": [{
            "ring_size": 6,
            "atom_indices": [0, 1, 2, 3, 4, 5],
            "dihedral_angles": [60.0, -60.0, 60.0, -60.0, 60.0, -60.0],
            "cremer_pople": {
                "Q": 0.7,
                "theta": 0.0,
                "phi": 0.0,
                "conformation": "chair"
            }
        }]
    }
    ```
    """
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol = Chem.AddHs(m.m)
        return _analyze_all_ring_conformations(mol, smiles, lang)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


def _analyze_all_ring_conformations(
    mol: Chem.Mol,
    smiles: str,
    lang: str = "ru",
) -> dict:
    """
    Генерирует 20 конформеров, вычисляет энергию каждого,
    определяет все уникальные конформации циклов с указанием
    наиболее стабильной (min energy) и наименее стабильной (max energy).

    Args:
        mol: RDKit молекула (с добавленными атомами водорода)
        smiles: исходная SMILES строка
        lang: язык ответа

    Returns:
        словарь с результатами анализа
    """
    # Находим все циклы (индексы атомов C/H в молекуле с водородами)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    if not rings:
        return {
            "smiles": smiles,
            "total_rings": 0,
            "rings": [],
            "message": get_text("no_rings_found", lang),
        }

    # Генерация 20 конформеров
    params = AllChem.ETKDGv3()
    params.numThreads = 0
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)

    if len(cids) == 0:
        raise ValueError(get_text("failed_to_generate_conformers", lang))

    # Оптимизация и вычисление энергии для каждого конформера
    conf_energies = {}
    for cid in cids:
        energy = 0.0
        try:
            energy = _optimize_and_energy(mol, cid)
        except Exception:
            try:
                # Fallback: UFF
                AllChem.UFFOptimizeMolecule(mol, confId=cid)
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                if ff:
                    energy = ff.CalcEnergy()
            except Exception:
                pass

        conf_energies[int(cid)] = round(float(energy), 4)

    # Анализируем каждый цикл для всех конформеров
    rings_data = []
    for ring in rings:
        ring_atoms = list(ring)
        unique_conformations = []
        seen_dihedrals = []  # Для проверки уникальности по диэдральным углам
        seen_energies = []  # Для проверки уникальности по энергии
        seen_pairs = []  # (форма, энергия)

        for cid in cids:
            if int(cid) not in conf_energies:
                continue

            conf_data = _get_ring_conformation_for_conf(mol, ring_atoms, int(cid))
            if conf_data and conf_data["cremer_pople"]:
                energy = conf_energies[int(cid)]
                
                # Проверяем уникальность по энергии (допуск 0.5 kcal/mol)
                is_unique_energy = True
                for e in seen_energies:
                    if abs(e - energy) < 0.5:
                        is_unique_energy = False
                        break

                # Проверяем уникальность по диэдральным углам (допуск 15 градусов)
                is_unique_dihedral = True
                for seen in seen_dihedrals:
                    if len(seen) == len(conf_data["dihedral_angles"]):
                        if all(abs(a - b) < 15.0 for a, b in zip(seen, conf_data["dihedral_angles"])):
                            is_unique_dihedral = False
                            break

                # Конформация уникальна, если уникальна по энергии ИЛИ по углам
                label = conf_data["cremer_pople"]["conformation"]
                is_new = not any(l == label and abs(e - energy) < 0.5 for l, e in seen_pairs)
                if is_new:
                    seen_pairs.append((label, energy))
                    seen_dihedrals.append(conf_data["dihedral_angles"])
                    seen_energies.append(energy)
                    conf_data["energy"] = energy
                    unique_conformations.append(conf_data)

        if unique_conformations:
            # Находим наиболее и наименее стабильные
            energies = [c["energy"] for c in unique_conformations]
            min_energy_idx = energies.index(min(energies))
            max_energy_idx = energies.index(max(energies))

            for i, c in enumerate(unique_conformations):
                c["is_most_stable"] = (i == min_energy_idx)
                c["is_least_stable"] = (i == max_energy_idx)

            rings_data.append({
                "ring_size": len(ring_atoms),
                "atom_indices": list(ring_atoms),
                "conformations": unique_conformations,
                "total_conformations_found": len(unique_conformations),
                "most_stable_conformation": unique_conformations[min_energy_idx]["cremer_pople"]["conformation"],
                "least_stable_conformation": unique_conformations[max_energy_idx]["cremer_pople"]["conformation"],
                "min_energy": min(energies),
                "max_energy": max(energies),
            })

    return {
        "smiles": smiles,
        "total_rings": len(rings_data),
        "rings": rings_data,
    }
