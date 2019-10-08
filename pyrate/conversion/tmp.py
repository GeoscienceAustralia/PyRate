str = "8 28 59.6906"
hrs, mins, sec = str.strip().split(" ")
time = int(hrs)*60*60 + int(mins)*60 + float(sec)
print(time)