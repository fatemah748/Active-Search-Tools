from rdkit import Chem
import sys

def readsmiles(smiles_path):
    with open(smiles_path) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
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
    with open("smiles_coreid-info.txt", "w+") as f:
        for key, value in dict.items():
            if value == None:
                continue
            else:
                f.write(str(key))
                f.write(", ")
                f.write(str(value))
                f.write("\n")

def main():
    # read in cores
    smiles_core_path = sys.argv[1]
    cores = readsmiles(smiles_core_path)

    # assign core ids to fragments
    coreids = {}
    n = 1
    for i in cores:
        coreids[n] = [i, 0]
        n += 1
        print(str(n - 1) + " " + i)
    # read in substituted smiles strings file
    subs_smile = sys.argv[2]
    subs_mol = readsmiles(subs_smile)
    # assign core ids to substituted strings
    subs_coreids = {}
    for s in subs_mol:
        m = Chem.MolFromSmiles(s)
        substituted_smiles = Chem.MolToSmiles(m, kekuleSmiles=True)
        for key, value in coreids.items():
            core_smi = Chem.CanonSmiles(value[0])
            p = Chem.MolFromSmiles(core_smi)
            if m.HasSubstructMatch(p):
                subs_coreids[substituted_smiles] = key
                coreids[key][1] += 1
                break
            else:
                subs_coreids[substituted_smiles] = None
        if subs_coreids[substituted_smiles] == None:
            print("No core found: " + str(Chem.MolToSmiles(m)))
    for key, value in coreids.items():
        if value[1] == 0:
            print("No deriviatives " + str(value[0]))
        else:
            print(str(value[0]) + ": " + str(value[1]))
    with open("info.txt", "w+") as f:
        for key, value in coreids.items():
            f.write(str(key) + " " + value[0] + " " + str(value[1]) + "\n")

    print(coreids)
    print(len(subs_coreids))
    #write_to_csv(subs_coreids)
# output strings with core id


if __name__== "__main__":
    main()
