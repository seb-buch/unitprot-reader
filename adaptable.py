
from collections import defaultdict, OrderedDict

class Entry(object):
    _defined_properties = {
        1: 'ID',
        2: 'sequence',
        3: 'name',
        4: 'source',
        5: 'Family',
        6: 'gene',
        7: 'stereo',
        8: 'N_terminus',
        9: 'C_terminus',
        10: 'PTM',
        11: 'cyclic',
        12: 'target',
        13: 'synthetic',
        14: 'antimicrobial',
        15: 'antibacterial',
        16: 'antigram_pos',
        17: 'antigram_neg',
        18: 'antifungal',
        19: 'antiyeast',
        20: 'antiviral',
        21: 'antiprotozoal',
        22: 'antiparasitic',
        23: 'antiplasmodial',
        24: 'antitrypanosomic',
        25: 'antileishmania',
        26: 'insecticidal',
        27: 'anticancer',
        28: 'antitumor',
        29: 'cell_line',
        30: 'tissue',
        31: 'cancer_type',
        32: 'anticancer_activity',
        33: 'anticancer_activity_test',
        34: 'antiangiogenic',
        35: 'toxic',
        36: 'cytotoxic',
        37: 'hemolytic',
        38: 'hemolytic_activity',
        39: 'hemolytic_activity_test',
        40: 'RBC_source',
        41: 'cell_cell',
        42: 'hormone',
        43: 'quorum_sensing',
        44: 'immunomodulant',
        45: 'antihypertensive',
        46: 'drug_delivery',
        47: 'cell_penetrating',
        48: 'tumor_homing',
        49: 'blood_brain',
        50: 'antioxidant',
        51: 'antiproliferative',
        52: 'DSSP',
        53: 'pdb',
        54: 'experim_structure',
        55: 'PMID',
        56: 'taxonomy',
        57: 'all_organisms',
        58: 'activity_viral',
        59: 'activity_viral_test',
        60: 'solubility',
        61: 'activity',
        62: 'activity_test',
        63: 'ribosomal',
        64: 'experimental',
        65: 'biofilm'
    }

    nproperties_used = 65

    def __init__(self, sequence, fasta_comment=None):
        self.sequence = sequence
        self.properties = defaultdict(list)
        self.nproperties = len(self.properties)
        self._properties_by_name = dict(zip(self._defined_properties.values(), self._defined_properties.keys()))
        if fasta_comment is not None:
            self.get_parameters_from_fasta(fasta_comment)

        self.properties[2] = [sequence]

    def to_fasta(self):
        fasta = ">"

        for i in range(self.nproperties_used):
            prop_value = self.properties[i+1]
            if len(prop_value) == 0:
                fasta += "_"
            else:
                fasta += ";".join(prop_value)

                # make sure the property ends with ";"
                fasta += ";"

            fasta += " "

        # remove trailing space
        fasta = fasta.strip(" ")

        fasta += "\n{}\n".format(self.sequence)

        return fasta

    def as_human_readable(self):
        content = "ADAPTABLE Entry:\n"

        for i in range(self.nproperties_used):
            prop_name = self._defined_properties[i+1]
            prop_value = self.properties[i+1]

            content += "  -> {}: {}\n".format(prop_name, "; ".join(prop_value))
        return content

    def get_parameters_from_fasta(self, line):
        if line[0] != ">":
            raise ValueError("FASTA line does not start with >")
        line = line[1:].split()

        for propid, value in enumerate(line):
            if value != "_":
                value = value.split(";")

                for val in value:
                    if val != "":
                        self.properties[propid+1].append(val)

    @property
    def name(self):
        try:
            name = self.properties[3][0]
        except IndexError:
            name = "UNKNOWN"

        return name

    def __repr__(self):
        return "< Adaptable entry: {} - sequence: {} >".format(self.name, self.sequence)

    def __str__(self):
        return "{}: {}".format(self.name, self.sequence)

    def __getitem__(self, item):
        if type(item) is int:
            return self.properties[item]
        try:
            return self.properties[self._properties_by_name[item]]
        except KeyError:
            raise KeyError("No such property: {}".format(item))


class Library(object):
    def __init__(self, fname=None, encoding="utf-8"):
        self.fname = fname
        self.encoding = encoding

        self.entries = OrderedDict()
        self.entries_list = []

    def read(self):
        character_replacements = [
            ("\\xa0", " "),  # non-break space
            ("\\x96", " "),  # start of guarded area
            ("\\xb5", "µ"),  # µ character
            ("\\xec", "µ"),  # µ character wrongly encoding as "ì"
            ("\\xb1", "±"),   # ± character
        ]
        if self.fname is None:
            raise ValueError("No filename defined. Please set the 'fname' attribute")

        with open(self.fname, encoding=self.encoding, errors="backslashreplace") as fp:
            sequence = None
            fasta_comment = None

            for lino, line in enumerate(fp):
                line = line.strip()
                if line == "":
                    continue

                # Try to replace html-oriented characters
                for bad, good in character_replacements:
                    line = line.replace(bad, good)

                if "\\" in line:
                    position = line.index("\\")
                    character = line[position:position+4]
                    print("Warning: the character '{}' is not encoded in utf-8 (representation in cp1252: '{}')".format(
                        character,
                        character.encode("cp1252").decode("unicode_escape")
                    ))

                if line[0] == ">":
                    fasta_comment = line
                elif fasta_comment is not None:
                    sequence = line
                    entry = Entry(sequence, fasta_comment)

                    fasta_comment = None

                    self.entries[sequence] = entry
                    self.entries_list.append(entry)

                    sequence = None

            print("{} lines read -> {} entries loaded\n".format(lino+1, len(self.entries)))

    def save(self, fname=None, verbose=True):
        if fname is None:
            fname = self.fname

        with open(fname, "w") as fp:
            for entry in self.entries.values():
                fp.write(entry.to_fasta())

        if verbose:
            print("Library saved to '{}'".format(fname))

    def __getitem__(self, item):
        if type(item) == int:
            return self.entries_list[item]
        else:
            return self.entries[item]