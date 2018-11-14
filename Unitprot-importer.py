#!/usr/bin/env python
import requests
from lxml import etree as ET
from lxml import objectify
from adaptable import Entry, Library
import sys
import os
import logging
from collections import defaultdict

# See https://www.uniprot.org/docs/keywlist for the complete list
BIOPROPERTIES_UNIPROT_KEYWORDS = {
    "antimicrobial": [
        "KW-0929",
        "KW-0081",
        "KW-0044",
    ],
    "antibacterial": [
        "KW-0081",
    ],
    "antifungal": [
        "KW-0295",
    ],
    "antiviral": [
        "KW-0930",
    ],
    "antitumor": [
        "KW-0043",
    ],
}


def get_uniprot_entries_from_query(query, verbose=True, max_length=50, reviewed=True):
    entries = []
    if verbose:
        print("Interrogating Unitprot with the following query: '{}' "
              "and a max sequence length of {} (reviewed={})... ".format(query,
                                                                         max_length, reviewed),
              end=""
              )
        sys.stdout.flush()

    response = ""
    cache_file = ".cache_uniprot_{}_{}_{}".format(query, max_length, reviewed)

    if os.path.isfile(cache_file):
        print("No need: cache file found and used!")
        with open(cache_file, "r") as fp:
            response = fp.read()
    else:
        try:
            if reviewed:
                reviewed = "yes"
            else:
                reviewed = "no"

            timeout = 60
            r = requests.get("https://www.uniprot.org/uniprot/?query={}+length:[1+TO+{}]+AND+reviewed:{}&columns=id&format=tab".format(query, max_length, reviewed),
                             timeout=timeout)
        except requests.exceptions.Timeout:
            print("ERROR! Unitprot did not respond within {} seconds!".format(timeout))
        else:
            if r.status_code != requests.codes.ok:
                if verbose:
                    print("Sorry Unitprot didn't like the query (Status code= {})".format(r.status_code))
            else:
                print("OK")
                response = r.text
                with open(cache_file, "w") as fp:
                    fp.write(response)

    for val in response.split('\n'):
        val = val.strip()

        if len(val) == 0:
            continue
        elif val == "Entry":
            continue

        entries.append(val)

    if verbose:
        print("Unitprot returned {} entries".format(len(entries)))
        sys.stdout.flush()
    return entries


def get_uniprot_entry_from_id(entry_id, verbose=False):
    content = None

    cache_dir = ".cache"
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    cache_file = os.path.join(cache_dir, "uniprot-{}".format(entry_id))

    if verbose:
        print("  Retrieving ID={} from Unitprot... ".format(entry_id), end="")

    if os.path.isfile(cache_file):
        with open(cache_file, 'r') as fp:
            content = ET.fromstring(fp.read().encode('utf-8'))[0]
        if verbose:
            print("No need: loading data from cache file")
    else:
        try:
            r = requests.get("https://www.uniprot.org/uniprot/{}.xml".format(entry_id),
                             timeout=5)
        except requests.exceptions.Timeout:
            print("ERROR! Unitprot does not respond!")
            return content

        if r.status_code != requests.codes.ok:
            if verbose:
                print("Sorry Unitprot didn't like the query (Status code= {}".format(r.status_code))
        else:
            response_text = r.text
            content = ET.fromstring(response_text.encode('utf-8'))[0]

            with open(cache_file, "w") as fp:
                fp.write(response_text)

            if verbose:
                print("OK")

    return content


def get_sequence_from_uniprot_xml(tree):
    """

    :param xml.etree.ElementTree.ElementTree tree: an unitprot entry formatted as xml
    :return: string
    """

    return tree.find("{*}sequence").text.replace("\n", "")


def get_tag(elem):
    return ET.QName(elem.tag).localname


def populate_entry_using_uniprot_xml(entry, tree, debug=False):
    KEYWORDS = {
        "antimicrobial":
            [
                "microbial",
            ],
        "antigram_pos":
            [
                "gram-positive",
            ],
        "antigram_neg":
            [
                "gram-negative",
            ],
        "antibacterial":
            [
                "bacterial",
            ],
        "antiviral":
            [
                "viral",
            ],
        "anticancer":
            [
                "cancer", "tumor", "anticancer", "antitumor",
            ],
        "antiprotozoal":
            [
                "protozoa",
            ],
        "antiplasmodial":
            [
                "plasmodi",
            ],
        "antiparasitic":
            [
                "parasit",
            ],
        "antitrypanosomic":
            [
                "trypanosom",
            ],
        "antileishmania":
            [
                "leishman",
            ],
        "insecticidal":
            [
                "insecticid",
            ],
        "toxic": [
            "toxic",
        ],
        "cytotoxic":
            [
                "cytotoxic",
            ],
        "antiangiogenic":
            [
                "antiangiogen",
            ],
        "hemolytic":
            [
                "hemolytic",
            ],
        "pdb":
            [
                "pdb",
            ],
        "PMID":
            [
                "PMID"
            ],
        "DSSP":
            [
                "DSSP"
            ]
    }

    potential_properties = defaultdict(list)

    entry_id = "UNKNOWN"

    def warning_cry(elem, element_type="element"):
        logger.warning("Uniprot Entry '{}' -> the following {} is ignored but could be interesting: {}".format(
            entry_id,
            element_type,
            ET.tounicode(elem)))
    # TODO: check space in text
    for elem in tree:
        tag = get_tag(elem)

        if tag == "accession":
            # get id
            entry["ID"].append("uniprot{}".format(elem.text))
            entry_id = elem.text
        elif tag == "gene":
            entry["gene"].append(elem[0].text)
        elif tag == "name":
            entry["name"].append(elem.text)

        elif tag == "protein":
            for subelem in elem[0]:
                entry["name"].append(subelem.text)
        elif tag == "organism":
            entry["source"].append(elem[0].text)
            entry["taxonomy"].append("NCBI:{}".format(elem[1].get("id"))) # TODO: use lineage instead
        elif tag == "dbReference":
            dbtype = elem.get("type")

            # Ignored databases
            if dbtype in [
                "GO",  # Gene Ontology
                "InterPro",
                "EC",  # ExPASy/Brenda
                "EMBL",  # European Nucleotide Archive
                "EnsemblBacteria",
                "OrthoDB",
                "Proteomes",
                "RefSeq",  # Reference genome sequences
                "PRINTS",  # Protein Motif fingerprint database
                "PATRIC",  # Pathosystems database
                "BioCyc",  # Pathway/Genome database
                "GeneID",
                "Gene3D",
                "PANTHER",  # Genome DB
                "SMART",
                "UniPathway",
                "KEGG",  # Genomic DB
                "HOGENOM",  # Genomic DB
                "OMA",  # Genome DB
                "UniGene",
                "MGI",  # Mouse Genome DB
                "UCSC",  # Genome Browser
                "Bgee",  # Gene expression DB
                "IntAct",  # Prot-Prot interaction DB
                "PeptideAtlas",  # Proteomic DB
                "PRIDE",  # Proteomic DB
                "FlyBase",  # Drosophila Genome DB
                "eggNOG",  # Genome DB
                "PDBsum",  # Ignored as redondant with PDB
                "PIR",  # Non structural Protein DB
                "STRING",  # Prot-Prot interaction
                "PaxDb",  # Protein Abundance Database
                "HOVERGEN",
                "TCDB",  # Transport Protein DB
                "MINT",  # Prot-Prot interaction
                "EvolutionaryTrace",  # Evolutionary DB
                "ArachnoServer",  # Prot DB dedicated to spider toxins
                "CAZy",  # Carbohydrate-active enzymes DB
                "Ensembl",  # Genome DB
                "PMAP-CutDB",  # Proteolytic DB
                "MaizeGDB",  # Maize Genome DB
                "iPTMnet",  # System biology DB
                "KO",  # Ortholog DB
                "Genevisible",
                "ExpressionAtlas",
                "SABIO-RK",  # Biochme kinetic DB
                "Allergome",  # Allergen DB
                "PRO",  # Protein ontology DB
                "DisProt",
                "InParanoid",  # Phylogenic DB
                "Araport",  # Arabido DB
                "ConoServer",  # Cone sanil toxin DB
                          ]:
                continue
            elif dbtype == "SUPFAM":  # Super Family of protein
                entry["ID"].append("supfam{}".format(elem.get("id")))
            elif dbtype == "Pfam":  # Family of protein
                entry["ID"].append("pfam{}".format(elem.get("id")))
            elif dbtype == "ProteinModelPortal":
                entry["ID"].append("pmp{}".format(elem.get("id")))
            elif dbtype == "TIGRFAMs":
                entry["ID"].append("tigrfams{}".format(elem.get("id")))
            elif dbtype == "HAMAP":
                entry["ID"].append("hamap{}".format(elem.get("id")))
            elif dbtype == "PROSITE":
                entry["ID"].append("prosite{}".format(elem.get("id")))
            elif dbtype == "PIRSF":
                entry["ID"].append("pirsf{}".format(elem.get("id")))
            elif dbtype == "CDD":
                entry["ID"].append("cdd{}".format(elem.get("id")))
            elif dbtype == "ProDom":
                entry["ID"].append("prodom{}".format(elem.get("id")))
            elif dbtype == "SMR":
                entry["ID"].append("smr{}".format(elem.get("id")))
            elif dbtype == "PDB":
                entry["pdb"].append(elem.get("id")) # TODO: also add it experment_structre
            else:
                warning_cry(elem, "Database ({})".format(dbtype))
        elif tag == "keyword":
            kwid = elem.get("id")

            for name, keywords in BIOPROPERTIES_UNIPROT_KEYWORDS.items():
                if kwid in keywords:
                    if len(entry[name]) == 0:
                        entry[name].append(name)
        else:
            if debug:
                print("DEBUG: Tag '{}' will be ignored".format(tag))

            flat_elem = ET.tounicode(elem).lower()
            for bioproperty, keywords in KEYWORDS.items():
                for keyword in keywords:
                    if keyword in flat_elem:
                        potential_properties[flat_elem].append(bioproperty)
                        break

    # Check potential properties
    log_entry = False
    for flat_elem, potential_properties in potential_properties.items():
        for bioproperty in potential_properties:
            if len(entry[bioproperty]) == 0:
                if bioproperty in ["DSSP", "pdb"]:
                    message = "structural information"
                elif bioproperty == "PMID":
                    message = "PMID"
                else:
                    message = "information about {} properties".format(bioproperty)

                logger.warning("Entry {} -> "
                               "Potential {} ignored but contained in the following element: {}".format(
                    entry_id,
                    message,
                    flat_elem
                ))
                log_entry = True

    if log_entry:
        logger.warning("Entry {} content:\n{}".format(
            entry_id,
            entry.as_human_readable()
        ))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process a Unitprot query and create a database accordingly.')
    parser.add_argument("query", help="Query to send to Unitprot server")
    parser.add_argument("--max-length", type=int, default=50,
                        help="Specify the maximum length (number of AA) for the Unitprot query")
    parser.add_argument("--verbose", action="store_true", help="Be verbose")
    parser.add_argument("--nonreviewed", action="store_false", help="Search for non-reviewed entries", dest="reviewed")
    parser.add_argument("--basename", help="Base name for the file where the database will be saved to", default="DATABASE")

    args = parser.parse_args()


    SILENT = not args.verbose

    logger = logging.getLogger('Uniport-importer')
    LOGFILE = 'uniprot_importer.log'
    os.remove(LOGFILE)
    hdlr = logging.FileHandler(LOGFILE)
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.WARNING)


    #print("Loading ADAPTABLE database:")

    #current_library = Library("../DATABASE")
    #current_library.read()

    unitprot_library = Library("{}_{}".format(args.basename, args.query))

    entries = get_uniprot_entries_from_query(args.query, verbose=args.verbose, max_length=args.max_length, reviewed=args.reviewed)

    counter_all = 0
    counter_new = 0
    errors = 0
    max_errors = 10
    try:
        for num, entry in enumerate(entries):
            print("\rProcessing entry {:5d}/{:5d}... ".format(num+1, len(entries)), end="")

            tree = get_uniprot_entry_from_id(entry)

            if tree is None:
                print("WARNING: Could not retrieve ID:{} from Uniprot"
                      ". It will be ignored".format(entry))
                errors += 1

                if errors > max_errors:
                    print("ERROR: Max errors ({}) reached... Aborting" .format(max_errors))
                    break
                continue

            counter_all += 1

            sequence = get_sequence_from_uniprot_xml(tree)
            ellipsed_sequence = sequence
            if len(ellipsed_sequence) > 50:
                ellipsed_sequence = ellipsed_sequence[:50] + "..."

            # try:
            #     raise KeyError
            #     entry = current_library.entries[sequence]
            #     if SILENT:
            #         print("UPDATE", end="")
            #     else:
            #         print("UPDATE entry:", end="")
            # except KeyError:
            #     entry = Entry(sequence)
            #
            #     if SILENT:
            #         print("NEW   ", end="")
            #     else:
            #         print("NEW entry:", end="")
            #     counter_new += 1

            entry = Entry(sequence)

            if SILENT:
                print("NEW   ", end="")
            else:
                print("NEW entry:", end="")
            counter_new += 1

            if not SILENT:
                print(" '{}'".format(ellipsed_sequence))

            populate_entry_using_uniprot_xml(entry, tree)

            unitprot_library.entries[sequence] = entry
    except KeyboardInterrupt:
        print("Keyboard interrupt triggered! Exiting")
    else:
        if SILENT:
            print("")
        unitprot_library.save()

    print("Summary: {} entries retrieved -> {} new entries (i.e. not already in ADAPTABLE)".format(counter_all,
                                                                                                   counter_new))
    print("All warning messages are saved in '{}'!".format(LOGFILE))






