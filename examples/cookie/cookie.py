def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(sys.argv) == 3:
        finname = sys.argv[1]
        foutname = sys.argv[2]
    hin = open(finname)
    hout = open(foutname, 'w')
    i = 0
    for line in hin:
        i = i + 1
        hout.write('Case #' + str(i) + ':\n')
        for row in solveMaze(line).mazeprintf():
            hout.write(row + '\n')
    hout.close()
    return 0


if __name__ == '__main__':
    main()
