# This Python file uses the following encoding: utf-8
import requests, re, time
import xml.etree.ElementTree as ET


class MagmaCalculator(object):
    def __init__(self):
        self.url = "https://magma.maths.usyd.edu.au/xml/calculator.xml"
        self.headers = { 'Content-Type': 'application/x-www-form-urlencoded' }

    def __parse_xml(self, response):
        root = ET.fromstring(response)
        lines_xml = root.iter("line")
        result = dict()
        k1 = 0; k2 = 0; lee = 0; bi_image = False
        val1 = 0; val2 = 0
        lines = [l.text for l in lines_xml]

        if len(lines) ==0:
            result['found'] = False
        else:
            for i in range(len(lines)):
                if i==0:
                    exponents = re.findall(r'\^(\d+)', lines[i])
                    k1 = int(exponents[0]); k2 = int(exponents[1])
                    result['k1'] = k1; result['k2'] = k2
                elif lines[i] is None:
                    continue
                elif lines[i].startswith('['):
                    continue
                elif lines[i].startswith('true'):
                    bi_image = True; result['bi_image'] = bi_image
                    match = re.search(r'\[(\d+),\s*(\d+)', lines[i])
                    if match:
                        val1 = int(match.group(1)); val2 = int(match.group(2))
                        result['bi_dim'] = val1; result['k1_k2_dim'] = val2
                    else:
                        print("thats bad there is no bi_dim etc...")
                elif lines[i].startswith('false'):
                        bi_image = False; result['bi_image'] = bi_image
                elif i == len(lines)-2:
                    lee = int(lines[i]); result['lee'] = lee
                else:
                    continue
            if bi_image:
                print(f"[ 4^{k1} 2^{k2}, {lee}] [{val1}, {val2}, {lee}]")
            else:
                print(f"[ 4^{k1} 2^{k2}, {lee}] [Not Linear]")
            result['found'] = True
        return result

    def _print_generator_matrix(self, G):
        d1, d2 = G.shape
        query_parts1 = ["Z4 := IntegerRing(4);", f"T := LinearCode<Z4, {d2} |"]
        rows_as_strings = []
        for row in G:
            row_str = "[" + ",".join(map(str, row)) + "]"
            rows_as_strings.append(row_str)
        query_body = ",\n".join(rows_as_strings)
        query_parts2 = ["T;", "HasLinearGrayMapImage(T);", "MinimumLeeWeight(T);"]
        query = "\n".join(query_parts1) + "\n" + query_body + ">;\n" + "\n".join(query_parts2)

        query = "input=" + query
        return query

    def calculate(self, G):

        query = self._print_generator_matrix(G)

        try:
            time.sleep(10)
            response2 = requests.post(self.url, data=query, headers=self.headers)
            if response2.status_code == 200:
                results = self.__parse_xml(response2.text) 
                results['query'] = query
            else:
                results = dict()
                results['found'] = False
                print(f"Sonuç Sayfası Hata kodu: {response2.status_code}")

            return results
        except Exception as e:
            print(f"Hata: {e}")


# if __name__ == "__main__":
#     pass
