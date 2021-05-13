"""
-------------------------------------------------------------------- 
---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ----------
--------------------------------------------------------------------

This file has been ported from the NRLMSISE-00 C source code package 
- release 20041227

The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
Doug Drob. Model is also available as a NRLMSISE-00 distribution 
package in FORTRAN (link to model FORTRAN: 
https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/)

Dominik Brodowski implemented and maintains this C version. You can
reach him at mail@brodo.de. 

Version Dated: 2019-07-09 1255 hrs

This is the Testing File for the NRLMSISE00 model implementation.
"""
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
    input = [nrlmsise_input() for _ in range(17)]
    flags = nrlmsise_flags()
    aph = ap_array()

    # Input Values
    for ind in range(0, 7):
        aph.a[ind] = 100

    flags.switches[0] = 0

    for ind in range(1, 24):
        flags.switches[ind] = 1

    for ind in range(0, 17):
        input[ind].doy = 172
        # without effect
        input[ind].year = 0
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

    return output, input


nrlmsise_output_values = {
    0: [
        666517.690495152,
        113880555.97522168,
        19982109.255734544,
        402276.3585712511,
        3557.464994515886,
        4.074713532757222e-15,
        34753.12399717142,
        4095913.2682930017,
        26672.73209335869,
        1250.5399435607994,
        1241.4161300191206,
    ],
    1: [
        3407293.223160913,
        158633336.9569168,
        13911173.65461115,
        326255.95095955464,
        1559.6181505012246,
        5.001845729072244e-15,
        48542.08463340254,
        4380966.712898625,
        6956.681955942268,
        1166.7543837572089,
        1161.7104518870424,
    ],
    2: [
        112376.72440379359,
        69341.3008676061,
        42.47105217477091,
        0.1322750141474932,
        2.6188484182321854e-05,
        2.7567723192688737e-18,
        20167.498543214315,
        5741.25593414718,
        23743.941519895943,
        1239.8921117166653,
        1239.890640133059,
    ],
    3: [
        54115543.79936674,
        191889344393.9309,
        6115825598224.634,
        1225201051740.1243,
        60232119730.84871,
        3.5844263041133343e-10,
        10598796.977405403,
        261573.6693705142,
        2.8198793559283335e-42,
        1027.3184649,
        206.8877764036055,
    ],
    4: [
        1851122.4861925277,
        147655483.7927462,
        15793562.282644961,
        263379.4977312314,
        1588.78139838393,
        4.809630239407452e-15,
        58161.667807874735,
        5478984.479068789,
        1264.4459417610085,
        1212.3961521212093,
        1208.1354252123917,
    ],
    5: [
        867309.5233906157,
        127886176.80141275,
        18225766.271717,
        292221.4190618247,
        2402.9624364237006,
        4.355865642644646e-15,
        36863.8924375055,
        3897275.503726964,
        26672.73209335869,
        1220.146417915032,
        1212.7120832118062,
    ],
    6: [
        577625.1216023244,
        69791386.93660195,
        12368135.598217027,
        249286.77154291022,
        1405.738674177843,
        2.4706513916631316e-15,
        52919.85567066641,
        1069814.1093665662,
        26672.732093358678,
        1116.3853760431516,
        1112.998568217311,
    ],
    7: [
        374030.41055076657,
        47827201.23611342,
        5240380.033324204,
        175987.46403906072,
        550.1648779569964,
        1.5718887392548444e-15,
        88967.75722935038,
        1979740.836232955,
        9121.814875991493,
        1031.247440714559,
        1024.848492213009,
    ],
    8: [
        674833.8766623624,
        124531526.04437314,
        23690095.410529852,
        491158.3154749823,
        4578.78109905442,
        4.564420245361171e-15,
        32445.947751610933,
        5370833.087086038,
        26672.73209335869,
        1306.0520420272921,
        1293.374040389534,
    ],
    9: [
        552860.0841645189,
        119804132.40413581,
        34957977.64558215,
        933961.8355028158,
        10962.547654934302,
        4.974543110322225e-15,
        26864.278562598087,
        4889974.232971402,
        28054.448371256636,
        1361.8680207849234,
        1347.3891837297017,
    ],
    10: [
        137548758418628.66,
        0,
        2.0496870442907566e19,
        5.498695433718813e18,
        2.4517331580283885e17,
        0.0012610656611185514,
        0,
        0,
        0,
        1027.3184649,
        281.4647576632156,
    ],
    11: [
        44274425876770.93,
        0,
        6.597567157737311e18,
        1.7699293414061885e18,
        7.891679955727485e16,
        0.0004059139375799179,
        0,
        0,
        0,
        1027.3184649,
        227.4179808272618,
    ],
    12: [
        2127828756207.1882,
        0,
        3.170790550354043e17,
        8.50627980943479e16,
        3792741116805988.5,
        1.9508222451756214e-05,
        0,
        0,
        0,
        1027.3184649,
        237.43891458772688,
    ],
    13: [
        141218354559.28513,
        0,
        2.104369643783158e16,
        5645392443377078.0,
        251714174941122.44,
        1.2947090159285671e-06,
        0,
        0,
        0,
        1027.3184649,
        279.555112954128,
    ],
    14: [
        12548844002.726698,
        0,
        1874532829219022.8,
        492305098078476.4,
        22396854138563.84,
        1.1476676715115365e-07,
        0,
        0,
        0,
        1027.3184649,
        219.07323136419572,
    ],
    15: [
        519647.740297288,
        127449407.29604626,
        48504498.69853357,
        1720837.9825749004,
        23544.86590544426,
        5.881940448651632e-15,
        25000.783910809296,
        6279209.8250188,
        26672.732093358678,
        1426.4116622824247,
        1408.607795553264,
    ],
    16: [
        42608597.487941206,
        124134201554.8743,
        4929561542488.143,
        1048406749092.832,
        49934650830.55505,
        2.9143035503087925e-10,
        8831228.592571594,
        225251.55086261546,
        2.4152459296489138e-42,
        1027.3184649,
        193.40710625766815,
    ],
}


@pytest.mark.parametrize("z", nrlmsise_output_values.keys())
def test_output(z):
    expected_output = [nrlmsise_output() for _ in range(17)]
    for ind in range(17):
        for memb in range(9):
            expected_output[ind].d[memb] = nrlmsise_output_values[z][memb]
        expected_output[ind].t[0] = nrlmsise_output_values[z][9]
        expected_output[ind].t[1] = nrlmsise_output_values[z][10]

    output_var, _ = test_gtd7()

    for element in range(9):
        assert output_var[z].d[element] == expected_output[z].d[element]
    assert output_var[z].t[0] == expected_output[z].t[0]
    assert output_var[z].t[1] == expected_output[z].t[1]


nrlmsise_input_output = {
    0: [172, 29000, 400, 60, -70, 16, 150, 150],
    1: [81, 29000, 400, 60, -70, 16, 150, 150],
    2: [172, 75000, 1000, 60, -70, 16, 150, 150],
    3: [172, 29000, 100, 60, -70, 16, 150, 150],
    4: [172, 29000, 400, 0, -70, 16, 150, 150],
    5: [172, 29000, 400, 60, 0, 16, 150, 150],
    6: [172, 29000, 400, 60, -70, 4, 150, 150],
    7: [172, 29000, 400, 60, -70, 16, 70, 150],
    8: [172, 29000, 400, 60, -70, 16, 150, 180],
    9: [172, 29000, 400, 60, -70, 16, 150, 150],
    10: [172, 29000, 0, 60, -70, 16, 150, 150],
    11: [172, 29000, 10, 60, -70, 16, 150, 150],
    12: [172, 29000, 30, 60, -70, 16, 150, 150],
    13: [172, 29000, 50, 60, -70, 16, 150, 150],
    14: [172, 29000, 70, 60, -70, 16, 150, 150],
    # 15: {},
    # 16: {},
}

@pytest.mark.parametrize("z", nrlmsise_input_output.keys())
def test_output(z):
    expected_input = [nrlmsise_input() for _ in range(15)]
    
    expected_input[z].doy = nrlmsise_input_output[z][0]
    expected_input[z].sec = nrlmsise_input_output[z][1]
    expected_input[z].alt = nrlmsise_input_output[z][2]
    expected_input[z].g_lat = nrlmsise_input_output[z][3]
    expected_input[z].g_long = nrlmsise_input_output[z][4]
    expected_input[z].lst = nrlmsise_input_output[z][5]
    expected_input[z].f107A = nrlmsise_input_output[z][6]
    expected_input[z].f107 = nrlmsise_input_output[z][7]

    _, input_var = test_gtd7()

    assert input_var[z].doy == expected_input[z].doy
    assert input_var[z].sec == expected_input[z].sec
    assert input_var[z].alt == expected_input[z].alt
    assert input_var[z].g_lat == expected_input[z].g_lat
    assert input_var[z].g_long == expected_input[z].g_long
    assert input_var[z].lst == expected_input[z].lst
    assert input_var[z].f107A == expected_input[z].f107A
    assert input_var[z].f107 == expected_input[z].f107  