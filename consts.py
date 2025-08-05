
# Known irreducible polynomials (binary):
# For full list until 10000 see https://shiftleft.com/mirrors/www.hpl.hp.com/techreports/98/HPL-98-135.pdf
IRREDUCIBLE_POLYS = {
    1: 0b11,  # x + 1
    2: 0b111,  # x^2 + x + 1
    3: 0b1011,  # x^3 + x + 1
    4: 0b10011,  # x^4 + x + 1
    5: 0b100101,  # x^5 + x^2 + 1
    6: 0b1000011,  # x^6 + x + 1
    7: 0b10000011,  # x^7 + x + 1
    8: 0b100011011,  # x^8 + x^4 + x^3 + x + 1
    9: 0b1000000101,  # x^9 + x^4 + 1
    10: 0b10000001001,  # x^10 + x^3 + 1
    11: 0b100000000101,  # x^11 + x^2 + 1
    12: 0b1000001010011,  # x^12 + x^6 + x^4 + x + 1
    13: 0b10000000011011,  # x^13 + x^4 + x^3 + x + 1
    14: 0b100010000000011,  # x^14 + x^10 + x^6 + x + 1
    15: 0b1000000000000011,  # x^15 + x + 1
}

QUBIT_NUM = 5

