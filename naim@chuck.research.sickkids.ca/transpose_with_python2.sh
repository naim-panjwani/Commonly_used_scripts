python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.strip().split() for l in sys.stdin.readlines() if l.strip()))))" < bowei.txt> bowei_transpose.txt
