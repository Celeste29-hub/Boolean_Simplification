import re
from itertools import compress, filterfalse
from collections import Counter, defaultdict

#Designed by Celeste Aurore Martin
class solve:
    def __init__(self, arg):
        #Define functions
        literals = lambda x: re.findall("\w'?", x)
        swap_dict = lambda x: dict(zip(x.values(), x.keys()))
        
        def combos(iterables):
            pool = literals(iterables)
            n = len(pool)
            for i in reversed(range(n)):
                out_1 = pool[i]
                yield out_1
                for j in range(i + 1, n):
                    out_1 = "".join(pool[i:j])
                    for k in reversed(range(j, n)):
                        yield out_1 + pool[k]
                
        atoms = sorted(Counter(re.findall("[a-z]", arg)).keys())
        if not atoms:
            raise Exception("None")
        self.atoms = atoms
        
        table = []
        atom_indices = defaultdict(str)
        for i in range(1 << len(atoms)):
            binary = bin(i)[2:].zfill(len(atoms))
            table += ["".join([x + "'" if not int(y) else x for x, y in zip(atoms, binary)])]
            atom_indices["1"] += "1"
            for key, value  in zip(atoms, binary):
                atom_indices[key] += value
        atom_indices = dict(zip(atom_indices.keys(), map(lambda x: int(x[::-1], 2), atom_indices.values())))
        
        #DNF minimization
        replacements = {
        r"\s|\'\'": "",
        r"([a-z'\)])([\(a-z])": r"\1&\2",
        r"([a-z])": lambda x: str(atom_indices[x.group(1)]),
        r"([0-9]+)\'": lambda x: str(int(x.group(1)) ^ atom_indices["1"]),
        r"(\w+)\&(\w+)": lambda x: str(int(x.group(1)) & int(x.group(2))),
        r"(\w+)\+(\w+)": lambda x: str(int(x.group(1)) | int(x.group(2))),
        r"\((\w+)\)": r"\1",
        }
        
        index = arg
        if re.findall("[^a-z\(\)\+\'\s]", arg):
            raise Exception("Check Argument")
        else:
            try:
                # Parses expression
                while any([re.findall(key, index) for key in replacements.keys()]):
                    for key, value in replacements.items():
                        while re.search(key, index):
                            index = re.sub(key, value, index)
                index = int(index)
                b_index = bin(int(index))[2:].zfill(1 << len(atoms))
                
                # Generate truth table
                table = dict(zip(table, map(int, b_index[::-1])))
                pdnf = list(compress(table.keys(), table.values()))
                self.table = table
                self.pdnf = pdnf
                
                # Specific cases
                if index == atom_indices["1"]:
                    self.solution = "1"
                elif index == 0:
                    self.solution = "0"
                elif index in atom_indices.values():
                    self.solution = swap_dict(atom_indices)[index]
                elif index ^ atom_indices["1"] in atom_indices.values():
                    self.solution = swap_dict(atom_indices)[index ^ atom_indices["1"]] + "'"
                elif [b_index[0], b_index[-1]] == ["0", "0"] and all([int(_) for _ in b_index[1:-1]]):
                    solution = [f"{''.join([atoms[i], atoms[(i + 1) % len(atoms)]])}'" for i in range(len(atoms))]
                    solution = ["".join(sorted(literals(_))) for _ in solution]
                    self.solution = " + ".join(sorted(solution))
                elif len(pdnf) == 1:
                    self.solution = pdnf[0]
                elif b_index.count("0") == 1:
                    solution = dict(zip(table.values(), table.keys()))[0]
                    solution = re.sub(r"(\w)", r"\1'", solution)
                    solution = re.sub("\'\'", "", solution)
                    solution = " + ".join(re.findall("\w'?", solution))
                    self.solution = solution
                else:
                    # Searches prime implicants
                    prime_implicants = defaultdict(list)
                    for key in pdnf:
                        for _ in combos(key):
                            prime_implicants[_] += [key]
                    
                    size = lambda x, y: len(y) < 1 << (len(atoms) - len(literals(x)))
                    prime_implicants = dict(filterfalse(lambda x: size(x[0], x[1]), prime_implicants.items()))
                    
                    groups = defaultdict(list)
                    for _ in prime_implicants.keys():
                        groups[len(literals(_))] += [_]
                    groups = {_: groups[_] for _ in sorted(groups)}
                    
                    # Find solution
                    dnf = []
                    temp = []
                    for _ in groups.values():
                        terms = dict(zip(_, map(lambda x: prime_implicants[x], _)))
                        while terms:
                            intersections  = lambda a, b: len(set(b)) - len(set(b) - set(a))
                            pref = sorted(terms.items(), key=lambda x: intersections(temp, x[1]))[0][0]
                            if set(terms[pref]) <= set(temp):
                                break
                            dnf += [pref]
                            temp += terms[pref]
                            del terms[pref]
                            if set(pdnf) == set(temp):
                                break
                        if set(pdnf) == set(temp):
                            break
                    
                    self.solution = " + ".join(sorted(sorted(dnf), key=lambda x: len(literals(x))))
            except Exception as e:
                raise Exception(f"{type(e).__name__}: {e}")

if __name__ == "__main__":
    while True:
            try:
                arg = solve(input("solve: "))
                print(f"F[{', '.join(arg.atoms)}] = {arg.solution}")
            except Exception as e:
                print(e)
            finally:
                input()
 
