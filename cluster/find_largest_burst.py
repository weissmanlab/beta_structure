import tskit
import argparse

parser = argparse.ArgumentParser(
                    prog='find_largest_burst')
parser.add_argument('--input', type=str, default="runs_structured/91")

args = parser.parse_args()

print("start")

ts = tskit.load(args.input)

print("tree sequence loaded")

### T OF LARGEST BURST
# For populations which have a coalescent event with arity > MIN_ARITY
# When is the event with largest arity?
trees_to_check = 20
min_arity = 2
highest_arity = 0
highest_arity_T = None
c=0
for tree in ts.trees():
    if c > trees_to_check:
        break
    for node in tree.nodes():
        if tree.num_children(node) >= min_arity and tree.num_children(node) > highest_arity:
            highest_arity = tree.num_children(node)
            highest_arity_T = tree.time(node)
    c+=1

with open(args.input + "_largest_burst", "w") as file:
    file.write(str(highest_arity) + "," + str(highest_arity_T) + "\n")

print("largest burst found")
