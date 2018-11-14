#!/usr/bin/env python
import requests
import xml.etree.ElementTree as ET

def get_entries_from_query(query, verbose=True):
    entries = []
    if verbose:
        print("Interrogating Unitprot with the following query: '{}'".format(query))
    try:
        r = requests.get("https://www.uniprot.org/uniprot/?query={}&columns=id&format=tab&limit=10".format(query),
                         timeout=5)
    except requests.exceptions.Timeout:
        print("ERROR! Unitprot does not respond!")

    if r.status_code != requests.codes.ok:
        if verbose:
            print("Sorry Unitprot didn't like the query (Status code= {}".format(r.status_code))
    else:
        for val in r.text.split('\n'):
            val = val.strip()

            if len(val) == 0:
                continue
            elif val == "Entry":
                continue

            entries.append(val)

        if verbose:
            print("Unitprot returned {} entries".format(len(entries)))
    return entries


def get_entry_from_id(entry_id, verbose=True):
    content = None
    if verbose:
        print("Retrieving ID={} from Unitprot... ".format(entry_id), end="")
    try:
        r = requests.get("https://www.uniprot.org/uniprot/{}.xml".format(entry_id),
                         timeout=5)
    except requests.exceptions.Timeout:
        print("ERROR! Unitprot does not respond!")

    if r.status_code != requests.codes.ok:
        if verbose:
            print("Sorry Unitprot didn't like the query (Status code= {}".format(r.status_code))
    else:
        content = ET.fromstring(r.text)

        if verbose:
            print("OK")
    return content


if __name__ == "__main__":
    entries = get_entries_from_query("antimicrobial")

    print(entries)

    for entry in entries:
        content = get_entry_from_id(entry)
        print(content)





