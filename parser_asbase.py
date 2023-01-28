import pandas as pd
import numpy as np
import os
from bs4 import BeautifulSoup
import h5py

string = h5py.string_dtype(encoding='utf-8')


def punctuation_pic_to_text(elem):
    pic_src = elem.img['src']
    if "None" in pic_src:
        return "Unknown"
    elif "True" in pic_src:
        return "Yes"
    elif "False" in pic_src:
        return "No"


def get_photostable(elem):
    return punctuation_pic_to_text(elem)


def get_thermostable(elem):
    return punctuation_pic_to_text(elem)


def get_solubility(elem):
    divs = elem.find_all("div", class_="col".split())
    res = []
    for div in divs:
        solution = div.get_text().strip()
        solubility = punctuation_pic_to_text(div)
        res.append(':'.join([solution, solubility]))
    return ",".join(res)


def parse_html(text):
    # pasrse html 2 dict
    data = {}
    soup = BeautifulSoup(text, features="lxml")
    card_groups = soup.find_all("div", {"class": "card-group"})
    skipped_headers = ["Chemical Structure"]

    for card in card_groups:
        elements = card.find_all("div", class_="card-header card-body".split())
        for ele in elements:
            class_name = ele['class'][0]
            if "header" in class_name:    # header
                header = ele.get_text(strip=True)

                if "Φ" in header:
                    header = header.replace("Φ", "phi ")
                if "(" in header:
                    header = header.replace("(", "")
                    header = header.replace(")", "")
                if "λ" in header:
                    header = header.replace("λ", "lambda ")
                if "τ" in header:
                    header = header.replace("τ", "tau ")
                if "ε" in header:
                    header = header.replace("ε", "epsilon ")
            elif "body" in class_name:    # body

                if "Photostable" in header:
                    bodystr = get_photostable(ele)
                elif "Thermostable" in header:
                    bodystr = get_thermostable(ele)
                elif "Solubility" in header:
                    bodystr = get_solubility(ele)
                else:
                    bodystr = ele.get_text().strip().replace("\n", "  ")
                if header not in skipped_headers:
                    data[header] = bodystr
    return data


def store_dat2h5(data, handler):
    id = data["ID"]
    group = handler.create_group(id)
    for key in data.keys():
        # dims = len(data[key])
        group.create_dataset(key, dtype=string, shape=(1))[...] = data[key]


handler = h5py.File("testdb.hdf5", "w")
files = os.listdir("./pages")
for text_file in files:
    text_file = "./pages/" + text_file

    text = open(text_file).read()

    data = parse_html(text)
    store_dat2h5(data, handler)
    print("Done on " + text_file)
handler.close()
