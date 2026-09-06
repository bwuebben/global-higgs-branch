"""Reproduce the degree-cutoff quantum magnetic alignment tables.

Run: sage -c "load('checks/check_x9_quantum_wall.sage')"
The GV coefficients are exact; these numerical roots and differences
between truncations are NOT bounds on the omitted instanton tail.
All covers of each retained GV term are included through Li2.
"""
import os
import mpmath as mp

previous_mirror_degree = os.environ.get('MIRROR_DEGREE')
try:
    os.environ['MIRROR_DEGREE'] = '8'
    load('checks/check_x9_mirror_periods.sage')
finally:
    if previous_mirror_degree is None:
        os.environ.pop('MIRROR_DEGREE',None)
    else:
        os.environ['MIRROR_DEGREE'] = previous_mirror_degree

def check_rounded(actual,expected):
    with mp.workdps(50):
        assert abs(mp.mpf(actual)-mp.mpf(expected))<mp.mpf('5e-13')

print('J*: label, polynomial root, degree-eight root, displacement')
for label,expected in [('F2','0.740632681385'),('FR','0.793233852071')]:
    result6 = truncated_wall(gv,6,label=label)
    result8 = truncated_wall(gv,8,label=label)
    check_rounded(result8[1],expected)
    with mp.workdps(50):
        assert abs(mp.mpf(result8[1])-mp.mpf(result6[1]))<mp.mpf('2e-20')
    print(label,*result8)

expected_approach = {
    QQ(1)/2:['0.772156594240','0.772156594879','0.772156594878'],
    QQ(1)/4:['0.763609038466','0.763611970047','0.763611987372'],
    QQ(1)/8:['0.756963017440','0.757213151227','0.757260907358'],
}
print('a=b=lambda=1, FR: epsilon, degree-four, degree-six, degree-eight roots')
for epsilon,expected in expected_approach.items():
    roots = [truncated_wall(gv,d,a=1,b=1,eps=epsilon,label='FR')[1]
             for d in [4,6,8]]
    for root,rounded in zip(roots,expected):
        check_rounded(root,rounded)
    print(epsilon,*roots)
print('ALL NUMERICAL TABLE CHECKS PASSED; cutoff agreement is not an error bound.')
