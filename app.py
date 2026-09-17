import urllib.parse
import urllib.request
from fastapi import FastAPI, HTTPException, Query, Depends
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from database import get_db, CachedName
from Molecule import Molecule

app = FastAPI(title="SMILES API Rdkit")

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

@app.get("/api/get_name")
def get_name(
    smiles: str = Query(..., description="SMILES to name"),
    db: Session = Depends(get_db)
):
    global PROXY_ADDRESS
    if not smiles.strip():
        raise HTTPException(status_code=400, detail="Empty SMILES")

    cached_record = db.query(CachedName).filter(CachedName.smiles == smiles).first()
    if cached_record:
        return {"smiles": smiles, "name": cached_record.name}

    safe_smiles = urllib.parse.quote(smiles, safe="")
    target_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_smiles}/iupac_name"

    try:
        opener = urllib.request.build_opener()
        if PROXY_ADDRESS.strip():
            proxy_host, proxy_port = PROXY_ADDRESS.split(":")
            import socks
            from sockshandler import SocksiPyHandler
            proxy_handler = SocksiPyHandler(
                socks.SOCKS5,
                proxy_host,
                int(proxy_port)
            )
            opener = urllib.request.build_opener(proxy_handler)

        req = urllib.request.Request(
            target_url,
            headers={'User-Agent': 'Mozilla/5.0'}
        )
        with opener.open(req, timeout=10) as resp:
            name = resp.read().decode('utf-8').strip()
            try:
                new_cache = CachedName(smiles=smiles, name=name)
                db.add(new_cache)
                db.commit()
            except IntegrityError:
                db.rollback() # Если параллельный поток успел записать раньше — откатываем и не падаем

            return {"smiles": smiles, "name": name}

    except urllib.error.HTTPError as e:
        db.rollback()
        if e.code == 404:
            raise HTTPException(status_code=404, detail="Name not found")
        raise HTTPException(status_code=e.code, detail=str(e))
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Proxy/Request error: {str(e)}")

@app.get("/api/render")
def render(
   smiles: str = Query(..., description="SMILES - of molecule"),
   draw_type: str = Query("just", description="Type: just, indexes, charges"),
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
       b64_str = Molecule.GetB64(mol_obj)
       return {
            "status": "ok",
            "smiles": smiles,
            "type": draw_type,
            "image_base64": b64_str,
       }
    except Exception as e:
       raise HTTPException(status_code=400, detail=str(e))
