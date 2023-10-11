using Polyhedra
using CDDLib

rep = hrep([HyperPlane([1, 1], 1)], [HalfSpace([0, -1], 0), HalfSpace([-1, 0], 0)])

halfspaces(rep)
hyperplanes(rep)

p = polyhedron(rep)

p2 = eliminate(p, [1])

hrep(p2)
rep

isempty(rep)

