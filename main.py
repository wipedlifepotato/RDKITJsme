from Molecule import Molecule
import sys , getopt
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import getopt, sys, uvicorn

host = '127.0.0.1'
port = 1235
PROXY_ADDRESS='127.0.0.1:9050'
def parse():
    global host, port
    global PROXY_ADDRESS
    def usage():
        print("""
-help -h
-port=1234 -p=1234
-host=127.0.0.1 -h
-s 127.0.0.1:9050
""")
    try:
        opts, args = getopt.getopt(sys.argv[1:], "up:h:s:", ["usage", "port=", "host=", "socks="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o == "-h":
            host = str(a)
        elif o in ('socks', 's'):
            PROXY_ADDRESS=str(a)
        elif o in ("-u", "--usage"):
            usage()
            sys.exit()
        elif o in ("-p", "--port"):
            port = int(port)
        else:
            assert False, "unhandled option"
    return True
import urllib.parse
import urllib.request
#TODO: to save a cache file, load in load. sigint signal
names = {}
def main():
    parse()
    app = FastAPI(title="SMILES API Rdkit")
    @app.get("/api/get_name")
    def get_name(
        smiles: str = Query(..., description="SMILES to name")
    ):
        global PROXY_ADDRESS,names
        if not smiles.strip():
            raise HTTPException(status_code=400, detail="Empty SMILES")
        safe_smiles = urllib.parse.quote(smiles, safe="")
        if smiles in names:
            return {"smiles": smiles, "name": names[smiles]}
        target_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_smiles}/iupac_name"
        try:
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
                names[smiles]=name
                return {"smiles": smiles, "name": name}
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise HTTPException(status_code=404, detail="Name not found")
            raise HTTPException(status_code=e.code, detail=str(e))
        except Exception as e:
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
    ###
    global host,port
    app.add_middleware(
    CORSMiddleware,
     allow_origins=["*"],
     allow_methods=["*"],
     allow_headers=["*"],
    )
    uvicorn.run(app, host=host, port=port)
if __name__ == "__main__":
    main()

