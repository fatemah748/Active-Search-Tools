from rdkit import Chem
import sys

def readsmiles(smiles_path):
    with open(smiles_path) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    content = [x.split(",") for x in content]
    return content

def fragsearch(m, p):
    matches = [set(x) for x in m.GetSubstructMatches(p)]
    frags = [set(y) for y in Chem.GetMolFrags(m)]
    for frag in frags:
            for match in matches:
                if match.issubset(frag):
                        return match
                else:
                    return False

def write_to_csv(dict):
    with open("iteration32-smiles-coreid-nm.txt", "w+") as f:
        for key, value in dict.items():
            f.write(str(key))
            f.write(", ")
            f.write(str(value[0]))
            f.write(",")
            f.write(str(value[1]))
            f.write("\n")

def main():
    # read in cores
    smiles_core_path = sys.argv[1]
    cores = readsmiles(smiles_core_path)

    # assign core ids to fragments
    coreids = {}
    n = 0
    for i in cores:
        n += 1
        coreids[n] = i
    print(coreids)
    # read in substituted smiles strings file
    subs_smile = sys.argv[2]
    subs_mol = readsmiles(subs_smile)
    print(len(subs_mol))
    # assign core ids to substitued strings
    subs_coreids = {}
    for s in subs_mol:
        m = Chem.MolFromSmiles(s[0])
        for key, value in coreids.items():
            p = Chem.MolFromSmiles(value[0])
            if m is not None:
                if m.HasSubstructMatch(p):
                    core_smiles = Chem.MolToSmiles(p)
                    substituted_smiles = Chem.MolToSmiles(m)
                    subs_coreids[substituted_smiles] = [key, s[1]]
    print(subs_coreids)

    write_to_csv(subs_coreids)


if __name__== "__main__":
    main()
