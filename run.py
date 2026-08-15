from CodeGenerator import CodeGenerator

generator = CodeGenerator()
thetas= ['a+2*b', 'b']
derivations=['2*b', '2*b']
grays= ['3*b', 'a+b']
# if x = 2 -> x^2 + x + 1 -->  n = 3
generator.set_parameters(u=1, x=7, z=4, thetas=thetas, derivations=derivations, grays=grays)
# below it wil run for thresold of 10.000 generated polynomials but if it finds 10 new or better codes will be ended before it hits the threshold. 
# new code means with this n,k1,k2 combination there is no code in http://quantumcodes.info/Z4/ list beside the found code
# better code means with this n,k1,k2 combination there is at least one code in http://quantumcodes.info/Z4/ list but the found code has the best min lee distance
# and also there is equal code which found code is equal to the best min lee distance code in http://quantumcodes.info/Z4/ list which written on different file.
generator.solve_for_random_polynomials(10, 100)


# for computation only one custom polynomial: 
# 6 & $3x^5+ (3u + 1)x^4+ x^3 + 3ux^2+ (3u + 2)x+ 3u + 1$* & $[12, 4^9 2^1, 2]$ & Not Linear \\
# coeffs = ['3','3*u+1','1','3*u','3*u+2','3*u+1']
# result_obj = generator.solve_for_custom_polynomial(coeffs)
# print(result_obj.get_debug_info())
# print(result_obj.condition)