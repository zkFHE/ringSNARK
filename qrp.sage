# Run me on https://sagecell.sagemath.org/
# Fill in the arrays v, w, y below with the 0/1 wiring predicates for your circuit, and this script will give you the coefficients for all v_k, w_k, y_k, as well as for t and h
q1 = 0xffffee001
q2 = 0xffffc4001
q3 = 0x1ffffe0001
q = q1*q2*q3
print("Modulus: ", q, "(", log(q, 2).n(), " bits )")

A = GF(q1) # We interpolate over exceptional set A = Z_q1
Ax = A['x']

def div_diff(xs, ys): # div_diff([(x_0, y_0), ..., (x_k, y_k)]) -> (div_diff([(x_1, y_1), ..., (x_k, y_k)] - div_diff([(x0, y0), ..., (x_{k-1}, y_{k-1})])) / (x_k - x_0)
    if len(xs) != len(ys):
        raise ValueError("xs and ys must have same length")
    if len(xs) == 0:
        raise ValueError("xs and ys must have non-zero length")
    if len(xs) == 1:
        return ys[0]
    res = div_diff(xs[1:], ys[1:]) - div_diff(xs[:-1], ys[:-1])
    return res / (xs[-1] - xs[0])

def interpolate(xs, ys):
    # cf. https://en.wikipedia.org/wiki/Newton_polynomial
    prod = 1
    res = div_diff(xs[0:1], ys[0:1]) # div_diff([(x_0, y_0)])
    for i in range(len(xs)-1):
        prod *= (x - xs[i]) # (x-x_i)
        res += div_diff(xs[:i+2], ys[:i+2]) * prod # div_diff([(x_0, y_0), ..., (x_{i+1}, y_{i+1})]) * prod
    return res


# For the circuit 
#      c6
#      |
#      *─────r6
#     / \c5
#    /   \
#   +     *──r5
#  / \   / \
# c1 c2 c3 c4

r5, r6 = var("r5 r6")
rs = [r5, r6]
v = [
    [0, 1], # v_1
    [0, 1], # v_2
    [1, 0],
    [0, 0], 
    [0, 0],
    [0, 0]
]
w = [
    [0, 0], 
    [0, 0], 
    [0, 0],
    [1, 0], 
    [0, 1],
    [0, 0]
]
y = [
    [0, 0], 
    [0, 0], 
    [0, 0],
    [0, 0], 
    [1, 0],
    [0, 1]
]

print("Coeffs: coeff_0,\tcoeff_1,\t...,\tcoeff_d")

v_poly = []
for i, vi in enumerate(v):
    poly = interpolate(rs, vi)
    v_poly.append(poly)
    print(f"v_{i+1}:\t", ",\t".join(map(str, poly.coefficients(x, sparse=False))))
print()
w_poly = []
for i, wi in enumerate(w):
    poly = interpolate(rs, wi)
    w_poly.append(poly)
    print(f"w_{i+1}:\t", ",\t".join(map(str, poly.coefficients(x, sparse=False))))
print()
y_poly = []
for i, yi in enumerate(y):
    poly = interpolate(rs, yi)
    y_poly.append(poly)
    print(f"y_{i+1}:\t", ",\t".join(map(str, poly.coefficients(x, sparse=False))))

print()
t = 1
for ri in rs:
    t *= (x-ri)
print("f:\t", ",\t".join(map(str, t.coefficients(x, sparse=False))))

V = 0
W = 0
Y = 0
cs = []
for i in range(len(v)):
    c_i = var(f"c{i+1}")
    cs.append(c_i)
    V += cs[i] * v_poly[i]
    W += cs[i] * w_poly[i]
    Y += cs[i] * y_poly[i]

numerator = V*W - Y

h, remainder = numerator.maxima_methods().divide(t)
print()
print("h:\t", ",\t".join(map(str, h.coefficients(x, sparse=False))))

print()
print("V:\t", ",\t".join(map(str, V.coefficients(x, sparse=False))))
print("W:\t", ",\t".join(map(str, W.coefficients(x, sparse=False))))
print("Y:\t", ",\t".join(map(str, Y.coefficients(x, sparse=False))))
# print("num:\t", ",\t".join(map(str, numerator.coefficients(x, sparse=False))))
