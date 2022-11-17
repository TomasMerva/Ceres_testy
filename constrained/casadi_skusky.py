from casadi import *


# Symbols/expressions
x = MX.sym('x')
y = MX.sym('y')
f = x**2 -x*y +y**2
g = -y+0.75

nlp = {}                 # NLP declaration
nlp['x']= vertcat(x,y)   # decision vars
nlp['f'] = f             # objective
nlp['g'] = g           # constraints

# Create solver instance
F = nlpsol('F','ipopt',nlp);

# Solve the problem using a guess
result = F(x0=[-1.8, 1.8],ubg=0,lbg=-float("inf"))
# result = F(x0=[-1.8, 1.8])

print(result['x'])