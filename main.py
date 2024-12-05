import re
from itertools import combinations, filterfalse
from collections import Counter, defaultdict

#Designed by Celeste Aurore Martin
class solve:
    def __init__(self, arg: str) -> str:
        if not isinstance(arg, str):
            raise Exception(f"{isinstance(arg, str) = }")
        arg = arg.lower()
        
        #Define functions
        literals = lambda x: re.findall("\w'?", x)
        swap_dict = lambda x: dict(zip(x.values(), x.keys()))
        
        def combos(iterables: str) -> 'generator':
            pool: list[str] = literals(iterables)
            n: int = len(pool)
            for i in range(1, n + 1):
                for _ in combinations(pool, i):
                    yield "".join(_)
                
        def table(string: str) -> 'generator':
            n: int = len(atoms)
            for i in reversed(range(1 << len(atoms))):
                binary: str = bin(i)[2:].zfill(n)
                term: str = ""
                for j in range(n):
                    term += atoms[j] 
                    term += "'" * int(binary[j])
                yield (term, int(string[i]))
                
        if not (atoms := sorted(Counter(re.findall("[a-z]", arg)).keys())):
            raise Exception("None")
        self.atoms: list[str] = atoms
        
        n: int = 1
        s: str = "1" * (1 << (len(atoms) - 1))
        atom_indices: dict = {}
        for i in range(len(atoms)):
            j = s[::-(1 << i)]
            j = j.zfill(len(j) << 1) * (1 << i)
            atom_indices[atoms[i]] = int(j[::-1], 2)
            n |= int(j[::-1], 2)
        atom_indices["1"] = n
        
        #DNF minimization
        replacements: dict = {
        r"\s|\'\'": "",
        r"([a-z'\)])([\(a-z])": r"\1&\2",
        r"([a-z])": lambda x: str(atom_indices[x.group(1)]),
        r"([0-9]+)\'": lambda x: str(int(x.group(1)) ^ atom_indices["1"]),
        r"(\w+)\&(\w+)": lambda x: str(int(x.group(1)) & int(x.group(2))),
        r"(\w+)\+(\w+)": lambda x: str(int(x.group(1)) | int(x.group(2))),
        r"\((\w+)\)": r"\1",
        }
        
        index: str = arg
        if re.findall("[^a-z\(\)\+\'\s]", arg):
            raise Exception("Check Argument")
        else:
            try:
                # Parses expression
                while any([re.findall(key, index) for key in replacements.keys()]):
                    for key, value in replacements.items():
                        while re.search(key, index):
                            index = re.sub(key, value, index)
                index: int = int(index)
                b_index: str = bin(int(index))[2:].zfill(1 << len(atoms))
                
                # Generate truth table
                pdnf: list[str] = list(key for key, value in table(b_index) if value)
                self.pdnf: str = " + ".join(pdnf)
                
                # Specific cases
                if index == atom_indices["1"]:
                    self.dnf: str = "1"
                elif index == 0:
                    self.dnf: str = "0"
                elif index in atom_indices.values():
                    self.dnf: str = swap_dict(atom_indices)[index]
                elif index ^ atom_indices["1"] in atom_indices.values():
                    self.dnf: str = swap_dict(atom_indices)[index ^ atom_indices["1"]] + "'"
                elif [b_index[0], b_index[-1]] == ["0", "0"] and all([int(_) for _ in b_index[1:-1]]):
                    solution: 'generator' = (f"{''.join([atoms[i], atoms[(i + 1) % len(atoms)]])}'" for i in range(len(atoms)))
                    solution: list[str] = ["".join(sorted(literals(_))) for _ in solution]
                    self.dnf: str = " + ".join(sorted(solution))
                elif len(pdnf) == 1:
                    self.dnf: str = pdnf[0]
                else:
                    # Searches prime implicants
                    prime_implicants: dict = defaultdict(list)
                    for key in pdnf:
                        for _ in combos(key):
                            prime_implicants[_] += [key]
                    
                    size = lambda x, y: len(y) < 1 << (len(atoms) - len(literals(x)))
                    prime_implicants = dict(filterfalse(lambda x: size(x[0], x[1]), prime_implicants.items()))
                    
                    groups: dict = defaultdict(list)
                    for _ in prime_implicants.keys():
                        groups[len(literals(_))] += [_]
                    groups = {_: groups[_] for _ in sorted(groups)}
                    
                    # Find solution
                    dnf: list[str] = []
                    temp: set[str] = set()
                    for _ in groups.values():
                        terms: dict = dict(zip(_, map(lambda x: prime_implicants[x], _)))
                        while terms:
                            intersections  = lambda a, b: len(set(b)) - len(set(b) - set(a))
                            pref: str = sorted(terms.items(), key=lambda x: intersections(temp, x[1]))[0][0]
                            if set(terms[pref]) <= temp:
                                break
                            dnf += [pref]
                            temp.update(terms[pref])
                            del terms[pref]
                            if set(pdnf) == temp:
                                break
                        if set(pdnf) == temp:
                            break
                    
                    self.dnf: str = " + ".join(sorted(sorted(dnf), key=lambda x: len(literals(x))))
            except Exception as e:
                raise Exception(f"{type(e).__name__}: {e}")

if __name__ == "__main__":
    while True:
            try:
                arg: str = solve(input("solve: "))
                print(f"F[{', '.join(arg.atoms)}] = {arg.dnf}")
            except Exception as e:
                print(e)
            finally:
                input()
 