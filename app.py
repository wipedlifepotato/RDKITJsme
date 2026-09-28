import os
import urllib.parse
import urllib.request
from fastapi import FastAPI, HTTPException, Query, Depends, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from database import get_db, CachedName
from Molecule import Molecule
from rdkit_router import router as rdkit_router, set_proxy as set_rdkit_proxy

from rdkit import Chem
from rdkit import DataStructs

from rdkit.Chem import AllChem, Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

app = FastAPI(
    title="SMILES API Rdkit",
    description="Хемоинформатический API на базе RDKit",
    version="2.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# Подключаем расширенный RDKit роутер
# include_router не работает — добавляем маршруты вручную
for route in rdkit_router.routes:
    app.routes.append(route)

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


@app.get("/api/similarity")
def calculate_similarity(
    smiles1: str = Query(..., description="First SMILES"),
    smiles2: str = Query(..., description="Second SMILES"),
):
    try:
        m1 = Molecule(smiles1).m
        m2 = Molecule(smiles2).m

        if m1 is None or m2 is None:
            raise ValueError("One or both SMILES are invalid")

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
def get_3d_sdf(smiles: str = Query(..., description="SMILES of molecule")):
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES")

        mol_3d = Chem.AddHs(m.m)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol_3d)

        sdf_block = Chem.MolToMolBlock(mol_3d)
        return Response(content=sdf_block, media_type="chemical/x-mdl-sdfile")
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/get_name")
def get_name(
    smiles: str = Query(..., description="SMILES to name"),
    db: Session = Depends(get_db),
):
    global PROXY_ADDRESS
    if not smiles.strip():
        raise HTTPException(status_code=400, detail="Empty SMILES")

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
                db.rollback()  # Если параллельный поток успел записать раньше — откатываем и не падаем

            return {"smiles": smiles, "name": name}

    except urllib.error.HTTPError as e:
        db.rollback()
        if e.code == 404:
            raise HTTPException(status_code=404, detail="Name not found")
        raise HTTPException(status_code=e.code, detail=str(e))
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Proxy/Request error: {str(e)}")


@app.get("/api/get_chiral")
def getChiralCenters(
    smiles: str = Query(..., description="SMILES - of molecule"),
):
    try:
        m = Molecule(smiles)
        if m.m is None:  # Небольшая проверка на валидность SMILES
            raise ValueError("Invalid SMILES string")

        # 2. Используем m.m вместо несуществующего mol
        centers = Chem.FindMolChiralCenters(m.m, includeUnassigned=True)

        return {"centers": centers, "centers_count": len(centers)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/get_properties")
def get_properties(smiles: str = Query(..., description="SMILES of molecule")):
    try:
        m = Molecule(smiles)
        if m.m is None:
            raise ValueError("Invalid SMILES string")

        mol = m.m
        return {
            "smiles": smiles,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),  # Коэффициент липофильности
            "hbd": Descriptors.NumHDonors(mol),  # Доноры водородных связей
            "hba": Descriptors.NumHAcceptors(mol),  # Акцепторы
            "tpsa": Descriptors.TPSA(mol),  # Полярная площадь поверхности
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/convert")
def convert_smiles(smiles: str = Query(..., description="SMILES to convert")):
    try:
        m = Molecule(smiles)
        if m.m is None:
            return {"valid": False, "error": "Invalid SMILES structure"}

        mol = m.m
        return {
            "valid": True,
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
            "isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "inchi": Chem.MolToInchi(mol),
            "inchikey": Chem.MolToInchiKey(mol),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/substructure_search")
def substructure_search(
    target_smiles: str = Query(..., description="Target molecule SMILES"),
    pattern_smarts: str = Query(..., description="SMARTS pattern to search for"),
):
    try:
        target = Molecule(target_smiles)
        pattern = Chem.MolFromSmarts(pattern_smarts)

        if target.m is None or pattern is None:
            raise ValueError("Invalid SMILES or SMARTS pattern")

        has_substruct = target.m.HasSubstructMatch(pattern)
        matches = target.m.GetSubstructMatches(pattern)  # Индексы совпавших атомов

        return {
            "has_match": has_substruct,
            "match_count": len(matches),
            "matched_atom_indices": matches,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/render")
def render(
    smiles: str = Query(..., description="SMILES - of molecule"),
    draw_type: str = Query("just", description="Type: just, indexes, charges, chiral"),
    highlight_chiral: bool = Query(False, description="Подсветить хиральные центры"),
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

        # Подсветка хиральных центров
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
