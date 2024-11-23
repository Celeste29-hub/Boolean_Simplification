import re
from itertools import combinations, compress, filterfalse
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
            for i in range(1, n + 1):
                for _ in combinations(pool, i):
                    yield "".join(_)
                
        atoms = sorted(Counter(re.findall("[a-z]", arg)).keys())
        if not atoms:
            raise Exception("None")
        self.atoms = atoms
        
        atom_indices = defaultdict(int)
        n = 1
        for i in range(len(atoms)):
            j = "1" * (1 << (len(atoms) - 1))
            k = j[::-(1 << i)]
            k = k.zfill(len(k) << 1) * (1 << i)
            atom_indices[atoms[i]] = int(k[::-1], 2)
            n |= int(k[::-1], 2)
        atom_indices["1"] = n
      
        table = []
        n = len(atoms)
        for i in reversed(range(1 << len(atoms))):
            binary = bin(i)[2:].zfill(n)
            term = ""
            for i in range(n):
                term += atoms[i] 
                term += "'" * int(binary[i])
            table += [term]
        
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
                self.pdnf = " + ".join(pdnf)
                
                # Specific cases
                if index == atom_indices["1"]:
                    self.dnf = "1"
                elif index == 0:
                    self.dnf = "0"
                elif index in atom_indices.values():
                    self.dnf = swap_dict(atom_indices)[index]
                elif index ^ atom_indices["1"] in atom_indices.values():
                    self.dnf = swap_dict(atom_indices)[index ^ atom_indices["1"]] + "'"
                elif [b_index[0], b_index[-1]] == ["0", "0"] and all([int(_) for _ in b_index[1:-1]]):
                    solution = [f"{''.join([atoms[i], atoms[(i + 1) % len(atoms)]])}'" for i in range(len(atoms))]
                    solution = ["".join(sorted(literals(_))) for _ in solution]
                    self.dnf = " + ".join(sorted(solution))
                elif len(pdnf) == 1:
                    self.dnf = pdnf[0]
                elif b_index.count("0") == 1:
                    solution = dict(zip(table.values(), table.keys()))[0]
                    solution = re.sub(r"(\w)", r"\1'", solution)
                    solution = re.sub("\'\'", "", solution)
                    solution = " + ".join(re.findall("\w'?", solution))
                    self.dnf = solution
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
                    
                    self.dnf = " + ".join(sorted(sorted(dnf), key=lambda x: len(literals(x))))
            except Exception as e:
                raise Exception(f"{type(e).__name__}: {e}")

if __name__ == "__main__":
    while True:
            try:
                arg = solve(input("solve: "))
                print(f"F[{', '.join(arg.atoms)}] = {arg.dnf}")
            except Exception as e:
                print(e)
            finally:
                input()
 
