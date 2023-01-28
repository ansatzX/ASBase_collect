from bs4 import BeautifulSoup
import requests
import re
import time

for i in range(1, 2435):
    url = "http://119.91.135.188:8080/sr/?srid=%d" % (i)
    r = requests.get(url)
    if r.status_code == 200:
        # save to disk
        open("pages/%d.html" % (i), "w").write(r.text)
        print("%d ok" % (i))
    else:
        print("%d fail" % (i))

    time.sleep(3)
