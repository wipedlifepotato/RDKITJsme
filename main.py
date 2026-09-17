import os
import sys
import getopt
import configparser
import uvicorn

from database import init_db
from app import set_proxy

def parse_config():
    host = '127.0.0.1'
    port = 1235
    proxy_address = '127.0.0.1:9050'
    database_url = 'sqlite:///cache.db'

    config_file = "config.ini"
    if os.path.exists(config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        
        if "server" in config:
            host = config["server"].get("host", host)
            port = config["server"].getint("port", port)
            proxy_address = config["server"].get("proxy", proxy_address)
            
        if "database" in config:
            database_url = config["database"].get("url", database_url)
    else:
        def usage():
            print("""
Usage:
  -h, --host=<host>        Bind host (default: 127.0.0.1)
  -p, --port=<port>        Bind port (default: 1235)
  -s, --socks=<proxy>      SOCKS proxy address (default: 127.0.0.1:9050)
  -d, --db=<url>           Database connection URL (default: sqlite:///cache.db)
  -u, --usage              Show this help message
            """)

        try:
            opts, args = getopt.getopt(
                sys.argv[1:], 
                "uh:p:s:d:", 
                ["usage", "host=", "port=", "socks=", "db="]
            )
        except getopt.GetoptError as err:
            print(err)
            usage()
            sys.exit(2)

        for o, a in opts:
            if o in ("-h", "--host"):
                host = str(a)
            elif o in ("-p", "--port"):
                port = int(a)
            elif o in ("-s", "--socks"):
                proxy_address = str(a)
            elif o in ("-d", "--db"):
                database_url = str(a)
            elif o in ("-u", "--usage"):
                usage()
                sys.exit(0)

    return host, port, proxy_address, database_url

def main():
    host, port, proxy_address, database_url = parse_config()
    
    init_db(database_url)
    set_proxy(proxy_address)
    
    uvicorn.run("app:app", host=host, port=port, reload=False)

if __name__ == "__main__":
    main()
