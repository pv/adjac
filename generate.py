import tempita
import sys


def main():
    process('adjac.f95.in', 'adjac.f95')


def process(src, dst):
    with open(src, 'rb') as f:
        text = f.read()
        tmpl = tempita.Template(text)
        out = tmpl.substitute()

    with open(dst, 'wb') as f:
        f.write(out)

    sys.exit(0)


if __name__ == "__main__":
    main()
