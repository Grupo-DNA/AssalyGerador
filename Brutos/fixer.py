import os

with open("gerar.txt", "r") as h1:
	r1 = h1.readlines()
	for l1 in r1:
		print(f"ID = {l1}")
		l1 = l1.replace("\n", "")
		with open("helper.txt", "r") as h2:
			r2 = h2.readlines()
			for l2 in r2:
				l2 = l2.replace("\n", "")
				with open(f"{l1}.txt", "r") as h3:
					with open(f"{l1}update.txt", "a") as h4:
						r3 = h3.readlines()
						for l3 in r3:
							sp=l3.split()
							if l2 == sp[0]:
								h4.write(f"{sp[0]}\t0\t0\t{sp[3]}{sp[4]}\n")
								#h4.write(f"{sp[0]}\t0\t0\t{sp[1]}\n")
		os.rename(f"{l1}.txt", f"{l1}up.txt")
		os.rename(f"{l1}update.txt", f"{l1}.txt")