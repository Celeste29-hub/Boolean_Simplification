#Designed by Celeste Aurore Martin 
import re
from itertools import combinations, filterfalse

#from collections import defaultdict
#For eductional purposes
class defaultdict(dict):
    def __init__(self, data_type) -> "datatype":
        self.data_type = data_type
    def __missing__(self, key) -> "datatype()":
        value = self.data_type()
        self[key] = value
        return value

class solve:
    def __init__(self, arg: str) -> str:
        if not isinstance(arg, str):
            raise Exception(f"{isinstance(arg, str) = }")
        self.arg = arg.lower()
        if not (atoms := sorted(set(re.findall("[a-z]", self.arg)))):
            raise Exception("None")
        if re.findall("[^a-z\(\)\+\'\s]", self.arg):
            raise Exception("Check Argument")
        
        self.atoms: list[str] = atoms
        self.index, self.atom_indices = self.parse()
        self.b_index: str = bin(int(self.index))[2:].zfill(1 << len(self.atoms))
        self.pdnf: set[str] = set(key for key, value in self.table(self.b_index) if value)
        if len(self.pdnf) == 1:
            self.dnf: str = self.pdnf[0]
        else:
            self.dnf: str = self.main()
        
    #Define functions
    literals = lambda self, x: re.findall("\w'?", x)
    swap_dict = lambda self, x: dict(zip(x.values(), x.keys()))
      
    #Generates a list of the combinations of a string of items
    def combos(self, iterables: list) -> iter:
        pool: list[str] = self.literals(iterables)
        n: int = len(pool)
        for i in range(1, n + 1):
            for _ in combinations(pool, i):
                yield "".join(_)
    
    #Generates a truth table for the given boolean expression
    def table(self, index: str) -> iter:
        n: int = len(self.atoms)
        for i in reversed(range(1 << len(self.atoms))):
            binary: str = bin(i)[2:].zfill(n)
            term: str = ""
            for j in range(n):
                term += self.atoms[j] 
                term += "'" * int(binary[j])
            yield (term, int(index[i]))
    
    #Parses the boolean expression to find the index number of the argument
    def parse(self) -> int:
        #Finds index values for the atoms used in the given boolean expression relative to our truth table
        n: int = 1
        s: str = "1" * (1 << (len(self.atoms) - 1))
        atom_indices: dict = {}
        for i in range(len(self.atoms)):
            j = s[::-(1 << i)]
            j = j.zfill(len(j) << 1) * (1 << i)
            atom_indices[self.atoms[i]] = int(j[::-1], 2)
            n |= int(j[::-1], 2)
        atom_indices["1"] = n
        
        replacements: dict = {
        r"\s|\'\'": "", #Removes white spaces and double quotations
        r"([a-z'\)])([\(a-z])": r"\1&\2", #ab == a&b
        r"([a-z])": lambda x: str(atom_indices[x.group(1)]), #Converts atoms to their index values relative to the given truth table
        r"([0-9]+)\'": lambda x: str(int(x.group(1)) ^ atom_indices["1"]), #Negation operation
        r"(\w+)\&(\w+)": lambda x: str(int(x.group(1)) & int(x.group(2))), #Conjunction operation
        r"(\w+)\+(\w+)": lambda x: str(int(x.group(1)) | int(x.group(2))), #Disjunction operation
        r"\((\w+)\)": r"\1", #Removes parentheses
        }
        
        #Parses expression
        index: str = self.arg
        while any([re.findall(key, index) for key in replacements.keys()]):
            for key, value in replacements.items():
                while re.search(key, index):
                    index = re.sub(key, value, index)
        return int(index), atom_indices
    
    def main(self) -> str:
        try:
            #Specific cases
            if self.index == self.atom_indices["1"]:
                return "1"
            elif self.index == 0:
                return "0"
            elif self.index in self.atom_indices.values():
                return self.swap_dict(self.atom_indices)[self.index]
            elif self.index ^ self.atom_indices["1"] in self.atom_indices.values():
                return self.swap_dict(self.atom_indices)[self.index ^ self.atom_indices["1"]] + "'"
            elif [self.b_index[0], self.b_index[-1]] == ["0", "0"] and all([int(_) for _ in self.b_index[1:-1]]):
                solution: iter = (f"{''.join([self.atoms[i], self.atoms[(i + 1) % len(self.atoms)]])}'" for i in range(len(self.atoms)))
                solution: list[str] = ["".join(sorted(self.literals(_))) for _ in solution]
                return " + ".join(sorted(solution))
            elif len(self.pdnf) == 1:
                return self.pdnf[0]
            
            size = lambda x, y: len(y) == 1 << (len(self.atoms) - len(self.literals(x)))
            intersections  = lambda a, b: len(set(b)) - len(set(b) - set(a))
            
            #Find solution
            def dnf() -> iter:
                #Searches prime implicants
                prime_implicants: dict = defaultdict(list)
                for key in self.pdnf:
                    for _ in self.combos(key):
                        prime_implicants[_] += [key]
                
                #Analyze literal size in implicants
                groups: dict = defaultdict(list)
                for key, value in prime_implicants.items():
                    if size(key, value):
                        groups[len(self.literals(key))] += [key]
                groups = {_: groups[_] for _ in sorted(groups)}
                
                #Tabulation Method
                temp: set[str] = set()
                for _ in groups.values():
                    terms: dict = {x: prime_implicants[x] for x in _}
                    while terms:
                        pref: str = sorted(terms.items(), key=lambda x: intersections(temp, x[1]))[0][0]
                        if set(terms[pref]) <= temp:
                            break
                        yield pref
                        temp.update(terms[pref])
                        del terms[pref]
                        if self.pdnf == temp:
                            return
            
            return " + ".join(sorted(sorted(dnf()), key=lambda x: len(self.literals(x))))
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
