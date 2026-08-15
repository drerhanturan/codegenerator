from datetime import datetime
from itertools import product
from random import choice, randint
import sys
from Distances import LeeDistance
from Magma import MagmaCalculator
import numpy as np
from Results import Result
from sympy import Matrix, Symbol, expand, lambdify, simplify, symbols, sympify


class CodeGenerator:

    def __init__(self):
        self.LeeDist = LeeDistance("lee_list.json") # collected from http://quantumcodes.info/Z4/ up to 2n = 40
        # Magma is a closed source licensed software but also give permission for a limited online calculation interface which used in this script 
        # to ask if a generated skew-cyclic code is linear or not beside that all computation is done locally by this class
        self.MagmaCalc = MagmaCalculator() # validation script for founded polynomials and codes from https://magma.maths.usyd.edu.au/calc/

    def set_parameters(self, u, x, z, thetas, derivations, grays, local=True):
        self.x_power = x + 1  # n = x + 1  
        self.u_power = u + 1
        self.u = symbols("u")
        self.z = z
        self.u_p_list = [self.u**k for k in range(self.u_power)]
        self.compute_locally = True

        coeff_symbols = [Symbol(chr(97 + i)) for i in range(self.u_power)]
        thetas_expr = self.simplify_coefficients(thetas)
        self._fast_theta_fns = [
            lambdify(coeff_symbols, t, modules="math") for t in thetas_expr
        ]
        derivations_expr = self.simplify_coefficients(derivations)
        self._fast_derive_fns = [
            lambdify(coeff_symbols, d, modules="math") for d in derivations_expr
        ]
        grays_expr = self.simplify_coefficients(grays)
        self._fast_gray_vec_fn = lambdify(
            coeff_symbols, grays_expr[: self.u_power], modules="numpy"
        )
        self.search_strategy = self.analyze_intersection_rules()

    def get_time(self, initial=None):
        if initial is None:
            initial = datetime.now()
            print("Started at: " + str(initial.strftime("%Y-%m-%d %H:%M:%S")))
            return initial
        else:
            now = datetime.now()
            passed = now - initial
            print("Ended at: " + str(now.strftime("%Y-%m-%d %H:%M:%S")))
            print("Time passed: " + str(passed.total_seconds()))
            return None

    # =========================================================================
    #  OVERLAPPING ANALYSIS
    # =========================================================================
    def _is_prime(self, n):
        if n < 2:
            return False
        for i in range(2, int(n**0.5) + 1):
            if n % i == 0:
                return False
        return True

    def analyze_intersection_rules(self):
        n = self.x_power
        gray_len = 2 * n

        rule1_octacode = gray_len % 8 == 0
        rule2_kerdock = gray_len % 4 == 0
        rule3_prime = self._is_prime(n)
        rule3_power_of_two = (n & (n - 1) == 0) and (n > 0)

   
        if rule1_octacode and rule3_power_of_two:
            mode = "SUPER_INTERSECTION_OCTACODE_RESONANCE"
            symmetric_only = True
            k2_even_only = True
            rank_tolerance = 1    
            aggressive_pruning = True
         
        elif rule1_octacode:
            mode = "OCTACODE_OPTIMAL_RANK_INTERSECTION"
            symmetric_only = True
            k2_even_only = True
            rank_tolerance = 1
            aggressive_pruning = True
        
        elif rule2_kerdock:
            mode = "KERDOCK_OPTIMAL_RANK"
            symmetric_only = False
            k2_even_only = True
            rank_tolerance = 1.5
            aggressive_pruning = False
         
        elif rule3_prime:
            mode = "PRIME_RESONANCE"
            symmetric_only = False
            k2_even_only = False
            rank_tolerance = 1
            aggressive_pruning = True
        else:
            mode = "STANDARD_FUNNEL"
            symmetric_only = False
            k2_even_only = False
            rank_tolerance = 2
            aggressive_pruning = False

        strategy = {
            "mode": mode,
            "rule1_octacode": rule1_octacode,
            "rule2_kerdock": rule2_kerdock,
            "rule3_prime": rule3_prime,
            "rule3_power_of_two": rule3_power_of_two,
            "symmetric_only": symmetric_only,
            "k2_even_only": k2_even_only,
            "rank_tolerance": rank_tolerance,
            "aggressive_pruning": aggressive_pruning,
            "optimal_target": n / 2.0,
        }

        print(
            f"\n[STRATEJİ ANALİZİ] x={self.x_power-1}, n={n}, Gray Uzunluğu={gray_len}"
        )
        print(f" Mod: {strategy['mode']}")
        print(
            f" Kurallar -> Octacode: {rule1_octacode}, Kerdock: {rule2_kerdock}, Asal: {rule3_prime}, 2^m: {rule3_power_of_two}\n"
        )
        return strategy

    # =========================================================================
    # Rule-based polynomial generator
    # =========================================================================
    def generate_algebraic_candidates(self):
 
        n = self.x_power
        units = [
            c for c in self.coefficients if self.get_ulist(c)[0] in [1, 3]
        ]
        all_coeffs = self.coefficients

         
        if self.search_strategy["symmetric_only"]:
            half_len = n // 2
            coeffs = [None] * n
 
            c0 = choice(units)
            coeffs[0] = c0 
            coeffs[n - 1] = self.reduce_polynomial(-self.get_theta(c0))

             
            for i in range(1, half_len):
                c_i = choice(all_coeffs)
                coeffs[i] = c_i
                coeffs[n - 1 - i] = self.reduce_polynomial(
                    -self.get_theta(c_i)
                )

            
            if n % 2 != 0:
                mid = half_len
                valid_mids = [
                    c
                    for c in all_coeffs
                    if self.reduce_polynomial(c + self.get_theta(c)) == 0
                ]
                coeffs[mid] = (
                    choice(valid_mids) if valid_mids else choice(all_coeffs)
                )

            return coeffs

        
        else:
            coeffs = [None] * n
            coeffs[0] = choice(units)
            coeffs[-1] = choice(units)
            for i in range(1, n - 1):
                coeffs[i] = choice(all_coeffs)
            return coeffs

    # =========================================================================
    # Rank check
    # =========================================================================
    def is_rank_promising(self, k1, k2):
        n = self.x_power
 
        if self.search_strategy["k2_even_only"] and (k2 % 2 != 0):
            return False

        rank_val = k1 + (k2 / 2.0)
        target = self.search_strategy["optimal_target"]
        tolerance = self.search_strategy["rank_tolerance"]

        if abs(rank_val - target) > tolerance:
            return False   

        return True

    def simplify_coefficients(self, expr_list):
        temp = [0] * len(expr_list)
        for i in range(len(expr_list)):
            temp[i] = simplify(expr_list[i]) % self.z
        return temp

    def generate_coefficients(self):
        if self.u_power < 2:
            raise NotImplementedError("u_power must be at least 2")
        u_powers = [self.u**i for i in range(self.u_power)]
        self.coefficients = []
        for coeffs in product(range(self.z), repeat=self.u_power):
            expr = 0
            for i in range(self.u_power):
                if coeffs[i] != 0:
                    expr += coeffs[i] * u_powers[i]
            if expr != 0:
                self.coefficients.append(self.reduce_polynomial(expr))

    def get_ulist(self, expr):
        ul = [0] * self.u_power
        expr = expand(expr)
        for term in expr.as_ordered_terms():
            coeff = term.as_coeff_exponent(self.u)[0]
            power = term.as_coeff_exponent(self.u)[1]
            ul[power] = int(coeff)
        return ul

    ##############################################################################
    #
    # This function is the custom function for the u^2 = 2, if you want to compute 
    # for another ring terms please consider to change the algorithm below for your 
    # custom quotient ring. Beside this function, for u exponents consider using 
    # u = 1 -> Z4 + uZ4
    # u = 2 -> Z4 + uZ4 + u^2Z4
    # u = 3 -> Z4 + uZ4 + u^2Z4 + u^3Z4 etc.
    # This is the only function that is not parameterized for skew-cyclic code searching
    # all other terms such as x, u, z and theta, derivation and gray.
    # for z it is tested for value 4, other ring values may not properly work because of 
    # custom rules of the generator.
    #
    ##############################################################################
    def reduce_polynomial(self, expr):
        expr = expand(simplify(expr))
        first_pass = 0
        for term in expr.as_ordered_terms():
            coeff = int(term.as_coeff_exponent(self.u)[0])
            power = term.as_coeff_exponent(self.u)[1]
            if power == 2:
                if coeff != 0:
                    coeff = 2 * coeff
                    power = 0
            first_pass += coeff * self.u**power
        first_pass = expand(simplify(first_pass))
        second_pass = 0
        for term in first_pass.as_ordered_terms():
            coeff = int(term.as_coeff_exponent(self.u)[0]) % self.z
            power = term.as_coeff_exponent(self.u)[1]
            second_pass += coeff * self.u**power
        return expand(second_pass)

    def multiply(self, p, q):
        return self.reduce_polynomial(p * q)

    def add(self, p, q):
        return self.reduce_polynomial(p + q)

    def get_theta(self, expr):
        ul = self.get_ulist(expr)
        res = 0
        for i, fn in enumerate(self._fast_theta_fns):
            res += (int(fn(*ul)) % self.z) * self.u_p_list[i]
        return self.reduce_polynomial(res)

    def get_derivation(self, expr):
        ul = self.get_ulist(expr)
        res = 0
        for i, fn in enumerate(self._fast_derive_fns):
            res += (int(fn(*ul)) % self.z) * self.u_p_list[i]
        return self.reduce_polynomial(res)

    def get_gray_map(self, n, M):
        M_array = np.array(
            [self.get_ulist(M[i, j]) for i in range(2 * n) for j in range(n)]
        )
        results = self._fast_gray_vec_fn(*M_array.T)
        stacked_results = np.stack(results, axis=-1)
        G_expanded = stacked_results.reshape(2 * n, 2 * n)
        G_expanded = G_expanded % self.z
        G_unique = np.unique(G_expanded, axis=0)
        G_non_zero_unique = G_unique[np.any(G_unique != 0, axis=1)]
        return G_non_zero_unique

    def get_generator_matrix(self, coefficients_list):
        c = len(coefficients_list)
        u_p_total = self.u_power
        base_rows = [None] * c
        base_rows[0] = coefficients_list[:]
        get_theta = self.get_theta
        get_derivation = self.get_derivation
        add = self.add
        multiply = self.multiply
        for k in range(1, c):
            prev = base_rows[k - 1]
            thetas = [get_theta(x) for x in prev]
            derivs = [get_derivation(x) for x in prev]
            shifted_thetas = [thetas[-1]] + thetas[:-1]
            base_rows[k] = [add(t, d) for t, d in zip(shifted_thetas, derivs)]
        u_factors = [self.u**p for p in range(u_p_total)]
        final_data = []
        for factor in u_factors:
            if factor == 1:
                final_data.extend(base_rows)
            else:
                for row in base_rows:
                    final_data.append([multiply(factor, val) for val in row])
        return Matrix(final_data)

    def print_generator_matrix(self, G):
        d1, d2 = G.shape
        query_parts1 = [
            "Z4 := IntegerRing(4);",
            f"T := LinearCode<Z4, {d2} |",
        ]
        rows_as_strings = []
        for row in G:
            row_str = "[" + ",".join(map(str, row)) + "]"
            rows_as_strings.append(row_str)
        query_body = ",\n".join(rows_as_strings)
        query_parts2 = [
            "T;",
            "HasLinearGrayMapImage(T);",
            "MinimumLeeWeight(T);",
        ]
        query = (
            "\n".join(query_parts1)
            + "\n"
            + query_body
            + ">;\n"
            + "\n".join(query_parts2)
        )

        query = "input=" + query
        return query

    def get_k1_k2_matrix(self):
        matrix = self.G.copy() % self.z
        rows, cols = matrix.shape
        pivot_row = 0

        for col in range(cols):
            if pivot_row >= rows:
                break
            odd_row = -1
            for r in range(pivot_row, rows):
                if matrix[r, col] % 2 == 1:
                    odd_row = r
                    break
            if odd_row != -1:
                matrix[[pivot_row, odd_row]] = matrix[[odd_row, pivot_row]]

                if matrix[pivot_row, col] == 3:
                    matrix[pivot_row] = (matrix[pivot_row] * 3) % self.z
                for r in range(rows):
                    if r != pivot_row and matrix[r, col] != 0:
                        factor = matrix[r, col]
                        matrix[r] = (
                            matrix[r] - factor * matrix[pivot_row]
                        ) % self.z
                pivot_row += 1
        k1 = pivot_row

        for col in range(cols):
            if pivot_row >= rows:
                break
            two_row = -1
            for r in range(pivot_row, rows):
                if matrix[r, col] == 2:
                    two_row = r
                    break
            if two_row != -1:
                matrix[[pivot_row, two_row]] = matrix[[two_row, pivot_row]]
                for r in range(rows):
                    if r != pivot_row and matrix[r, col] == 2:
                        matrix[r] = (matrix[r] - matrix[pivot_row]) % self.z
                pivot_row += 1
        k2 = pivot_row - k1
        return k1, k2, matrix[:pivot_row]

    def z4_vector_lee_weight(self, z4_vectors):
        a = z4_vectors % 2
        b = (z4_vectors // 2) % 2
        bin_img = np.hstack([b % 2, (a + b) % 2])
        return np.sum(bin_img, axis=1)

    def calculate(self, k1, k2, G_reduced, target_best_lee, max_samples=100_000):
        num_basis = G_reduced.shape[0]
        if num_basis == 0:
            return {"found": False}

        total_cw_count = (4**k1) * (2**k2)
        min_lee = float("inf")

        if total_cw_count <= max_samples:
            scalars_k1 = [list(range(4)) for _ in range(k1)]
            scalars_k2 = [[0, 1] for _ in range(k2)]
            all_combos = np.array(
                list(product(*(scalars_k1 + scalars_k2))), dtype=np.int8
            )
            z4_cw = (all_combos @ G_reduced) % self.z
            weights = self.z4_vector_lee_weight(z4_cw)
            non_zeros = weights[weights > 0]
            min_lee = int(np.min(non_zeros)) if len(non_zeros) > 0 else 0

        else:
            low_weight_combos = []
            for i in range(num_basis):
                for val in range(1, 4) if i < k1 else [1]:
                    vec = np.zeros(num_basis, dtype=np.int8)
                    vec[i] = val
                    low_weight_combos.append(vec)
                for j in range(i + 1, min(i + 5, num_basis)):
                    vec = np.zeros(num_basis, dtype=np.int8)
                    vec[i] = 1
                    vec[j] = 1
                    low_weight_combos.append(vec)
            combos_arr = np.array(low_weight_combos, dtype=np.int8)
            z4_cw = (combos_arr @ G_reduced) % self.z
            weights = self.z4_vector_lee_weight(z4_cw)
            non_zeros = weights[weights > 0]
            if len(non_zeros) > 0:
                min_lee = int(np.min(non_zeros))

            if not min_lee <= target_best_lee:
                rand_k1 = (
                    np.random.randint(
                        0, 4, size=(max_samples, k1), dtype=np.int8
                    )
                    if k1 > 0
                    else np.empty((max_samples, 0), dtype=np.int8)
                )
                rand_k2 = (
                    np.random.randint(
                        0, 2, size=(max_samples, k2), dtype=np.int8
                    )
                    if k2 > 0
                    else np.empty((max_samples, 0), dtype=np.int8)
                )
                rand_combos = np.hstack([rand_k1, rand_k2])
                z4_rand_cw = (rand_combos @ G_reduced) % self.z
                rand_weights = self.z4_vector_lee_weight(z4_rand_cw)
                rand_non_zeros = rand_weights[rand_weights > 0]
                if len(rand_non_zeros) > 0:
                    sampled_min = int(np.min(rand_non_zeros))
                    if sampled_min < min_lee:
                        min_lee = sampled_min

        return {
            "found": True,
            "k1": k1,
            "k2": k2,
            "lee": min_lee,
            "total_codewords": total_cw_count,
            "is_exact": total_cw_count <= max_samples,
        }

    def solve_for_custom_polynomial(self, coefficients_list):
        clean_coeffs = []
        for c in coefficients_list:
            if isinstance(c, str):
                clean_coeffs.append(sympify(c, locals={"u": self.u}))
            else:
                clean_coeffs.append(c)
        coefficients_list = clean_coeffs

        n = len(coefficients_list)
        M = self.get_generator_matrix(coefficients_list)
        self.G = self.get_gray_map(n, M)
        self.G = self.G[np.any(self.G != 0, axis=1)]
        self.G = np.unique(self.G, axis=0)
        self.k, self.m = self.G.shape

       
        k1, k2, G_reduced = self.get_k1_k2_matrix()

        
        if not self.is_rank_promising(k1, k2):
            return Result(
                coefficients_list,
                {"found": False, "reason": "rank_outside_optimal_band", "query": self.print_generator_matrix(G_reduced)},
            )

        best_lee = self.LeeDist.get_min_lee_dist_for_k1_k2(n * 2, k1, k2)

        
        candidate = self.calculate(k1, k2, G_reduced, best_lee)

        if candidate["found"]:
            result = dict()
            
            if candidate["lee"] > best_lee:
                print(self.print_generator_matrix(G_reduced))
                result = self.MagmaCalc.calculate(G_reduced)
                if result["found"]:
                    lee = result["lee"]
                    if lee > best_lee:
                        if best_lee == 0:
                            result["condition"] = "new"
                        else:
                            result["condition"] = "better"
                result["n"] = self.x_power * self.u_power
                result_obj = Result(coefficients_list, result)
                return result_obj
            elif candidate["lee"] == best_lee:
                result["found"] = True
                result["k1"] = k1
                result["k2"] = k2
                result["lee"] = candidate["lee"]
                result["n"] = self.x_power * self.u_power
                result["bi_image"] = False
                result["condition"] = "equal"
                result["query"] = self.print_generator_matrix(G_reduced)
                result_obj = Result(coefficients_list, result)
                return result_obj
            else:
                result["found"] = True
                result["k1"] = k1
                result["k2"] = k2
                result["lee"] = candidate["lee"]
                result["n"] = self.x_power * self.u_power
                result["bi_image"] = False
                result["condition"] = "below"
                result["query"] = self.print_generator_matrix(G_reduced)
                result_obj = Result(coefficients_list, result)
                return result_obj

        else:
            return Result(coefficients_list, {"found": False, "query": self.print_generator_matrix(G_reduced)})

    def solve_for_random_polynomials(self, sample_size, threshold):
        useful = open("usefuls.txt", "w", encoding="utf-8")
        equal = open("equals.txt", "w", encoding="utf-8")

        self.generate_coefficients()

        try:
            useful_count = 0
            total_attempts = 0
            initial = self.get_time()

            while useful_count < sample_size:
                total_attempts += 1

                print(
                    f"\r[SEARCHING] Try: {total_attempts}/{threshold} | Founded usefuls: {useful_count}/{sample_size} | Search Strategy: {self.search_strategy['mode']}",
                    end="",
                    flush=True,
                )

                if total_attempts > threshold:
                    print(
                        f"\nTotal attempts exceeded the threshold. Stopping at: {total_attempts} attempts."
                    )
                    self.get_time(initial)
                    break

                
                coeffs = self.generate_algebraic_candidates() 
                result_obj = self.solve_for_custom_polynomial(coeffs)

                if result_obj.found:
                    if result_obj.condition in ["new", "better"]:
                        useful_count += 1
                        useful.write(result_obj.get_debug_info())
                        useful.write(result_obj.condition + "\n\n\n")
 
                        print(
                            f"\n\n>>> Founded a new or better code! ({result_obj.condition.upper()})"
                        )
                        print(result_obj.get_debug_info())
                        print("-" * 50)

                    elif result_obj.condition == "equal":
                        equal.write(result_obj.get_debug_info())
                        equal.write(result_obj.condition + "\n\n\n")

        finally:
            print("\n[COMPLETED] Files are closed. Please check for the results.")
            useful.close()
            equal.close()