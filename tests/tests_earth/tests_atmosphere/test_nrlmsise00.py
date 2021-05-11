import astropy
import numpy as np
import pytest

from poliastro.earth.atmosphere.nrlmsise00 import (
    ap_array,
    gtd7,
    nrlmsise_flags,
    nrlmsise_input,
    nrlmsise_output,
)


def test_gtd7():
    output = [nrlmsise_output() for _ in range(17)]
    # output = nrlmsise_output()
    input = [nrlmsise_input() for _ in range(17)]
    # input = nrlmsise_input()
    flags = nrlmsise_flags()
    aph = ap_array()

    # input Values
    for ind in range(0, 7):
        aph.a[ind] = 100

    flags.switches[0] = 0

    for ind in range(1, 24):
        flags.switches[ind] = 1

    for ind in range(0, 17):
        # input[ind] = nrlmsise_input(0, 172, 29000, 400, 60, -70, 16, 150, 150, 4)
        input[ind].doy = 172
        input[ind].year = 0
        # without effect
        input[ind].sec = 29000
        input[ind].alt = 400
        input[ind].g_lat = 60
        input[ind].g_long = -70
        input[ind].lst = 16
        input[ind].f107A = 150
        input[ind].f107 = 150
        input[ind].ap = 4

    input[1].doy = 81
    input[2].sec = 75000
    input[2].alt = 1000
    input[3].alt = 100
    input[10].alt = 0
    input[11].alt = 10
    input[12].alt = 30
    input[13].alt = 50
    input[14].alt = 70
    input[16].alt = 100
    input[4].g_lat = 0
    input[5].g_long = 0
    input[6].lst = 4
    input[7].f107A = 70
    input[8].f107 = 180
    input[9].ap = 40

    input[15].ap_a = aph
    input[16].ap_a = aph

    # Evaluate 0 to 14
    for ind in range(0, 15):
        gtd7(input[ind], flags, output[ind])

    # Evaluate 15 and 16
    flags.switches[9] = -1
    for ind in range(15, 17):
        gtd7(input[ind], flags, output[ind])

    # # Output Type 1
    # for ind in range(0, 17):
    #     print("\n")
    #     for ind2 in range(0, 9):
    #         print(f"{output[ind].d[ind2]}")
    #     print(f"{output[ind].t[0]}")
    #     print(f"{output[ind].t[1]} \n")
    # # /* DL omitted */

    # # /* output type 2 */
    # for i in range(0, 3):
    #     print("\n")
    #     print("\nDAY   ")
    #     for j in range(0, 5):
    #         print(f"{input[i*5+j].doy}")
    #     print("\nUT    ")
    #     for j in range(0, 5):
    #         print(f"{input[i*5+j].sec}")
    #     print("\nALT   ")
    #     for j in range(0, 5):
    #         print(f"{input[i * 5 + j].alt}")
    #     print("\nLAT   ")
    #     for j in range(0, 5):
    #         print(f"{input[i * 5 + j].g_lat}")
    #     print("\nLONG  ")
    #     for j in range(0, 5):
    #         print(f"{input[i * 5 + j].g_long}")
    #     print("\nLST   ")
    #     for j in range(0, 5):
    #         print(f"{input[i * 5 + j].lst}")
    #     print("\nF107A ")
    #     for j in range(0, 5):
    #         print(f"{input[i * 5 + j].f107A}")
    #     print("\nF107  ")
    #     for j in range(0, 5):
    #         print(f"{input[i * 5 + j].f107}")
    #     print("\n\n")
    #     print("\nTINF  ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].t[0]}")
    #     print("\nTG    ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].t[1]}")
    #     print("\nHE    ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[0]}")
    #     print("\nO     ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[1]}")
    #     print("\nN2    ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[2]}")
    #     print("\nO2    ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[3]}")
    #     print("\nAR    ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[4]}")
    #     print("\nH     ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[6]}")
    #     print("\nN     ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[7]}")
    #     print("\nANM 0 ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[8]}")
    #     print("\nRHO   ")
    #     for j in range(0, 5):
    #         print(f"{output[i * 5 + j].d[5]}")
    #     print("\n")
    # print("\n")

    # /* output type 1 */
    for i in range(17):
        print("\n", end="")
        for j in range(9):
            print("%E " % output[i].d[j], end="")
        print("%E " % output[i].t[0], end="")
        print("%E " % output[i].t[1])
        # /* DL omitted */

    # /* output type 2 */
    for i in range(3):
        print("\n", end="")
        print("\nDAY   ", end="")
        for j in range(5):
            print("         %3i" % input[i * 5 + j].doy, end="")
        print("\nUT    ", end="")
        for j in range(5):
            print("       %5.0f" % input[i * 5 + j].sec, end="")
        print("\nALT   ", end="")
        for j in range(5):
            print("        %4.0f" % input[i * 5 + j].alt, end="")
        print("\nLAT   ", end="")
        for j in range(5):
            print("         %3.0f" % input[i * 5 + j].g_lat, end="")
        print("\nLONG  ", end="")
        for j in range(5):
            print("         %3.0f" % input[i * 5 + j].g_long, end="")
        print("\nLST   ", end="")
        for j in range(5):
            print("       %5.0f" % input[i * 5 + j].lst, end="")
        print("\nF107A ", end="")
        for j in range(5):
            print("         %3.0f" % input[i * 5 + j].f107A, end="")
        print("\nF107  ", end="")
        for j in range(5):
            print("         %3.0f" % input[i * 5 + j].f107, end="")

        print("\n\n", end="")

        print("\nTINF  ", end="")
        for j in range(5):
            print("     %7.2f" % output[i * 5 + j].t[0], end="")
        print("\nTG    ", end="")
        for j in range(5):
            print("     %7.2f" % output[i * 5 + j].t[1], end="")
        print("\nHE    ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[0], end="")
        print("\nO     ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[1], end="")
        print("\nN2    ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[2], end="")
        print("\nO2    ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[3], end="")
        print("\nAR    ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[4], end="")
        print("\nH     ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[6], end="")
        print("\nN     ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[7], end="")
        print("\nANM   ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[8], end="")
        print("\nRHO   ", end="")
        for j in range(5):
            print("   %1.3e" % output[i * 5 + j].d[5], end="")
        print("\n")


if __name__ == "__main__":
    test_gtd7()
