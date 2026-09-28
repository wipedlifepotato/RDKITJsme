import os
import json
import urllib.parse
import urllib.request
import urllib.error
from fastapi import FastAPI, HTTPException, Query, Depends, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from database import get_db, CachedName, NameToSmilesCache
from Molecule import Molecule
from rdkit_router import router as rdkit_router, set_proxy as set_rdkit_proxy
from i18n import get_text, get_lang

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import inchi  # ИСПРАВЛЕНИЕ: Отдельный импорт для генерации InChI

app = FastAPI(
    title="SMILES API Rdkit",
    description="Хемоинформатический API на базе RDKit",
    version="2.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# Инициализация БД (создание таблиц при первом запуске)
from database import init_db
init_db("sqlite:///./chem.db")

# Подключаем расширенный RDKit роутер
# include_router не работает с текущей версией FastAPI — добавляем маршруты вручную
for route in rdkit_router.routes:
    app.routes.append(route)

# ИСПРАВЛЕНИЕ: Создаем директорию, чтобы FastAPI не падал при старте, если ее нет
os.makedirs("StaticFiles", exist_ok=True)

# Раздаём статические файлы из StaticFiles (JSME и др.)
app.mount(
    "/StaticFiles",
    StaticFiles(directory="StaticFiles"),
    name="StaticFiles",
)


@app.get("/", response_class=HTMLResponse, include_in_schema=False)
async def serve_index():
    """Отдаёт основной JSME-клиент по корневому пути."""
    index_path = os.path.join(
        os.path.dirname(__file__), "StaticFiles", "index_jsme.html"
    )
    if os.path.exists(index_path):
        with open(index_path, "r", encoding="utf-8") as f:
            return f.read()
    return "<h1>SMILES API Rdkit</h1><p>index_jsme.html не найден</p>"


app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

PROXY_ADDRESS = "127.0.0.1:9050"


def set_proxy(proxy: str):
    global PROXY_ADDRESS
    PROXY_ADDRESS = proxy
    set_rdkit_proxy(proxy)


@app.get("/api/similarity")
def calculate_similarity(
    smiles1: str = Query(..., description=get_text("param_smiles1")),
    smiles2: str = Query(..., description=get_text("param_smiles2")),
    lang: str = Depends(get_lang),
):
    try:
        m1 = Molecule(smiles1).m
        m2 = Molecule(smiles2).m

        if m1 is None or m2 is None:
            raise ValueError(get_text("invalid_smiles", lang))

        fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, 2, nBits=2048)

        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return {
            "smiles1": smiles1,
            "smiles2": smiles2,
            "tanimoto_similarity": round(similarity, 4),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/get_3d_sdf")
def get_3d_sdf(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol_3d = Chem.AddHs(m.m)
        res = AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        if res == -1:
            raise ValueError(get_text("failed_to_embed_molecule", lang))
        AllChem.MMFFOptimizeMolecule(mol_3d)

        sdf_block = Chem.MolToMolBlock(mol_3d)
        return Response(content=sdf_block, media_type="chemical/x-mdl-sdfile")
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/get_name")
def get_name(
    smiles: str = Query(..., description=get_text("param_smiles_to_name")),
    lang: str = Depends(get_lang),
    db: Session = Depends(get_db),
):
    if not smiles.strip():
        raise HTTPException(status_code=400, detail=get_text("empty_smiles", lang))

    cached_record = db.query(CachedName).filter(CachedName.smiles == smiles).first()
    if cached_record:
        return {"smiles": smiles, "name": cached_record.name}

    safe_smiles = urllib.parse.quote(smiles, safe="")
    target_url = (
        f"https://cactus.nci.nih.gov/chemical/structure/{safe_smiles}/iupac_name"
    )

    try:
        opener = urllib.request.build_opener()
        if PROXY_ADDRESS.strip():
            proxy_host, proxy_port = PROXY_ADDRESS.split(":")
            import socks
            from sockshandler import SocksiPyHandler

            proxy_handler = SocksiPyHandler(socks.SOCKS5, proxy_host, int(proxy_port))
            opener = urllib.request.build_opener(proxy_handler)

        req = urllib.request.Request(target_url, headers={"User-Agent": "Mozilla/5.0"})
        with opener.open(req, timeout=10) as resp:
            name = resp.read().decode("utf-8").strip()
            try:
                new_cache = CachedName(smiles=smiles, name=name)
                db.add(new_cache)
                db.commit()
            except IntegrityError:
                db.rollback()  # Если параллельный поток успел записать раньше — откатываем

            return {"smiles": smiles, "name": name}

    except urllib.error.HTTPError as e:
        db.rollback()
        if e.code == 404:
            raise HTTPException(status_code=404, detail=get_text("name_not_found", lang))
        raise HTTPException(status_code=e.code, detail=str(e))
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"{get_text('proxy_request_error', lang)}: {str(e)}")


@app.get("/api/get_chiral")
def getChiralCenters(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        centers = Chem.FindMolChiralCenters(m.m, includeUnassigned=True)

        return {"centers": centers, "centers_count": len(centers)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/get_properties")
def get_properties(
    smiles: str = Query(..., description=get_text("param_smiles")),
    lang: str = Depends(get_lang),
):
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError(get_text("invalid_smiles", lang))

        mol = m.m

        # Точная масса через pyOpenMS
        exact_mass = None
        try:
            import pyopenms as oms
            formula = rdMolDescriptors.CalcMolFormula(mol)
            ef = oms.EmpiricalFormula(formula)
            exact_mass = ef.getMonoWeight()
        except Exception:
            pass

        return {
            "smiles": smiles,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": Descriptors.MolWt(mol),
            "exact_mass": exact_mass,
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/convert")
def convert_smiles(
    smiles: str = Query(..., description=get_text("param_smiles_to_convert")),
    lang: str = Depends(get_lang),
):
    try:
        m = Molecule(smiles)
        if m.m is None:
            return {"valid": False, "error": get_text("invalid_smiles_structure", lang)}

        mol = m.m
        return {
            "valid": True,
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
            "isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "inchi": inchi.MolToInchi(mol),        # ИСПРАВЛЕНИЕ: Использование подмодуля inchi
            "inchikey": inchi.MolToInchiKey(mol),  # ИСПРАВЛЕНИЕ: Использование подмодуля inchi
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/substructure_search")
def substructure_search(
    target_smiles: str = Query(..., description=get_text("param_target_smiles")),
    pattern_smarts: str = Query(..., description=get_text("param_pattern_smarts")),
    lang: str = Depends(get_lang),
):
    try:
        target = Molecule(target_smiles)
        pattern = Chem.MolFromSmarts(pattern_smarts)

        if target.m is None or pattern is None:
            raise ValueError(get_text("invalid_smiles_or_smarts", lang))

        has_substruct = target.m.HasSubstructMatch(pattern)
        matches = target.m.GetSubstructMatches(pattern)

        return {
            "has_match": has_substruct,
            "match_count": len(matches),
            "matched_atom_indices": matches,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/render")
def render(
    smiles: str = Query(..., description=get_text("param_smiles")),
    draw_type: str = Query("just", description=get_text("param_draw_type")),
    highlight_chiral: bool = Query(False, description=get_text("param_highlight_chiral")),
    lang: str = Depends(get_lang),
):
    try:
        m = Molecule(smiles)
        type_map = {
            "just": Molecule.DrawType.JUST,
            "indexes": Molecule.DrawType.INDEXES,
            "charges": Molecule.DrawType.CHARGES,
        }
        selected_type = type_map.get(draw_type.lower(), Molecule.DrawType.JUST)
        mol_obj = m.Get(selected_type)

        highlight_atoms = []
        if highlight_chiral or draw_type.lower() == "chiral":
            chiral_centers = Chem.FindMolChiralCenters(m.m, includeUnassigned=True)
            highlight_atoms = [c[0] for c in chiral_centers]

        if highlight_atoms:
            img = Draw.MolToImage(
                mol_obj,
                size=(500, 500),
                highlightAtoms=highlight_atoms,
                highlightAtomColors={a: (1.0, 0.2, 0.2) for a in highlight_atoms},
            )
            from io import BytesIO
            import base64

            buffer = BytesIO()
            img.save(buffer, format="PNG")
            b64_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        else:
            b64_str = Molecule.GetB64(mol_obj)

        return {
            "status": "ok",
            "smiles": smiles,
            "type": draw_type,
            "image_base64": b64_str,
            "chiral_centers": highlight_atoms
            if highlight_chiral or draw_type.lower() == "chiral"
            else [],
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/name_to_smiles")
def name_to_smiles(
    name: str = Query(..., description="Название молекулы (IUPAC или тривиальное)"),
    lang: str = Depends(get_lang),
    db: Session = Depends(get_db),
):
    """
    Конвертирует название молекулы (IUPAC или тривиальное) в SMILES.

    Использует OPSIN API (https://opsin.ch.cam.ac.uk) для конвертации
    названий в SMILES. При неудаче использует Cactus (NCI/NIH).
    Результаты кешируются в SQLite для ускорения повторных запросов.

    **Параметры:**
    - `name` — название молекулы (например, "Aspirin", "2-acetoxybenzoic acid")
    - `lang` — язык ответа (ru, en, uk, es)

    **Возвращает:**
    - `name` — исходное название
    - `smiles` — полученный SMILES
    - `source` — источник результата ("opsin", "cactus" или "cache")

    **Пример запроса:**
    ```
    GET /api/name_to_smiles?name=Aspirin&lang=ru
    ```

    **Пример ответа:**
    ```json
    {
        "name": "Aspirin",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "source": "opsin"
    }
    ```

    **Примечание:**
    При повторном запросе того же названия результат будет взят из кеша.
    """
    if not name.strip():
        raise HTTPException(status_code=400, detail=get_text("empty_name", lang))

    # Проверяем кеш
    cached = db.query(NameToSmilesCache).filter(NameToSmilesCache.name == name).first()
    if cached:
        return {"name": name, "smiles": cached.smiles, "source": "cache"}

    # Пробуем OPSIN API
    try:
        safe_name = urllib.parse.quote(name, safe="")
        target_url = f"https://opsin.ch.cam.ac.uk/opsin/{safe_name}.json"

        opener = urllib.request.build_opener()
        if PROXY_ADDRESS.strip():
            proxy_host, proxy_port = PROXY_ADDRESS.split(":")
            import socks
            from sockshandler import SocksiPyHandler

            proxy_handler = SocksiPyHandler(socks.SOCKS5, proxy_host, int(proxy_port))
            opener = urllib.request.build_opener(proxy_handler)

        req = urllib.request.Request(target_url, headers={"User-Agent": "Mozilla/5.0"})
        with opener.open(req, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            smiles = data.get("smiles", "")

        if smiles:
            # Сохраняем в кеш
            try:
                new_cache = NameToSmilesCache(name=name, smiles=smiles)
                db.add(new_cache)
                db.commit()
            except IntegrityError:
                db.rollback()

            return {"name": name, "smiles": smiles, "source": "opsin"}
    except Exception:
        pass

    # Fallback на Cactus
    try:
        safe_name = urllib.parse.quote(name, safe="")
        target_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_name}/smiles"

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

        if smiles:
            try:
                new_cache = NameToSmilesCache(name=name, smiles=smiles)
                db.add(new_cache)
                db.commit()
            except IntegrityError:
                db.rollback()

            return {"name": name, "smiles": smiles, "source": "cactus"}
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"{get_text('proxy_request_error', lang)}: {str(e)}",
        )

    raise HTTPException(status_code=404, detail=get_text("name_not_found", lang))
