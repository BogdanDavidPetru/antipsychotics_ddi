
nr_list = {12867, 12710, 12401, 13256, 13273, 13552, 13665, 13784, 13791, 13841, 13523, 4872, 4888, 8927, 11736, 12093, 12273, 12518, 12958, 13213, 13382, 13403, 13554, 13557, 13676, 372, 11540, 13420, 6077, 9266, 14185, 14651, 16021, 5766, 17056, 12543}

without_interaction = {'DB13665', 'DB00509', 'DB05766', 'DB13557', 'DB13523', 'DB09226', 'DB06077', 'DB13382', 'DB13256', 'DB13784', 'DB04888', 'DB11540', 'DB00372', 'DB12518', 'DB13403', 'DB16021', 'DB13552', 'DB12543', 'DB12093', 'DB12401', 'DB17056', 'DB12273', 'DB11376', 'DB13841', 'DB12693', 'DB14651', 'DB13420', 'DB12710', 'DB13273', 'DB08927', 'DB12867', 'DB12958', 'DB09223', 'DB14185', 'DB04872', 'DB13554', 'DB13213', 'DB13676', 'DB13791'}

found = 0

not_found = []

for nr in nr_list:
    s = str(nr)
    s_found = False
    for drug in without_interaction:
        if drug.__contains__(s):
            found = found + 1
            s_found = True
    if not s_found:
        not_found.append(s)

print(found)
print(not_found)

with_info_but_no_interaction =[]

for drug in without_interaction:
    found = False
    for nr in nr_list:
        s = str(nr)
        if drug.__contains__(s):
            found = True
            break
    if not found:
        with_info_but_no_interaction.append(drug)

print(with_info_but_no_interaction)
