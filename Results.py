# This Python file uses the following encoding: utf-8

class Result:
    def __init__(self, coeffs, result_dict):
        self.coeff_list = coeffs
        self.query = result_dict['query'][6:]
        self.found = result_dict['found']
        if self.found:
            self.condition = result_dict['condition']
            self.n = result_dict['n']
            self.k1 = result_dict['k1']
            self.k2 = result_dict['k2']
            self.lee = result_dict['lee']
            self.is_linear = result_dict['bi_image']
            if self.is_linear:
                self.bi_dim = result_dict['bi_dim']
                self.k1_k2_dim = result_dict['k1_k2_dim']


    def get_debug_info(self):
        if self.found:
            debug_text = self.get_latex_format()+"\n"
            debug_text = debug_text + self.query + "\n"
            return debug_text
        else:
            return "No result found..."

    def get_latex_format(self):
        if self.found:
            order = len(self.coeff_list)
            polynomial = "$"
            for i in range(len(self.coeff_list)):
                order = order - 1

                cur_coeff = str(self.coeff_list[i]).replace('*', '')

                if '+' in cur_coeff and order > 0:
                    cur_coeff = "(" + cur_coeff + ")"

                if order > 9:
                    x_base = "{" + str(order) + "}"
                else:
                    x_base = str(order)

                if order > 1:
                    if cur_coeff == "1":
                        polynomial = polynomial + "x^" + x_base + " + "
                    else:
                        polynomial = polynomial + cur_coeff + "x^" + x_base + " + "
                elif order == 1:
                    if cur_coeff == "1":
                        polynomial = polynomial + "x" + x_base + " + "
                    else:
                        polynomial = polynomial + cur_coeff + "x" + " + "
                else:
                    polynomial = polynomial + cur_coeff


            polynomial =polynomial + f"$ & $[{self.n}, "
            if self.k1>9:
                polynomial = polynomial + "4^{"+str(self.k1)+"}, "
            else:
                polynomial = polynomial + f"4^{self.k1}, "

            if self.k2>9:
                polynomial = polynomial + "2^{"+str(self.k2)+"}, "+str(self.lee)+"]$ & "
            else:
                polynomial = polynomial + f"2^{self.k2}, "+str(self.lee)+"]$ & "

            if self.is_linear:
                polynomial = polynomial + f"$[{self.bi_dim}, {self.k1_k2_dim}, {self.lee}]$ \\\\"
            else:
                polynomial = polynomial + "Not Linear \\\\"

            return polynomial
        else:
            return "No result found..."





# if __name__ == "__main__":
#     pass
