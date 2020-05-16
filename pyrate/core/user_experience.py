import pathlib


def delete_tsincr_files(params):

    out_dir = pathlib.Path(params["outdir"])
    for file_path in out_dir.iterdir():
        if "tsincr" in str(file_path):
            file_path.unlink()


def factorise_integer(n, memo={}, left=2):
    if (n, left) in memo:
        return memo[(n, left)]
    if left == 1:
        return n, [n]
    i = 2
    best = n
    bestTuple = [n]
    while i * i <= n:
        if n % i == 0:
            rem = factorise_integer(n / i, memo, left - 1)
            if rem[0] + i < best:
                best = rem[0] + i
                bestTuple = [i] + rem[1]
        i += 1

    # handle edge case when only one processor is available
    if bestTuple == [4]:
        return 2, 2

    if len(bestTuple) == 1:
        bestTuple.append(1)

    return int(bestTuple[0]), int(bestTuple[1])
