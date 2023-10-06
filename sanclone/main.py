# -*- coding: utf-8 -*-
import click

WELCOME = """
Welcome to San Clone 👋  a molecular cloning agent 🧬.
Give it an instruction like "Clone NADH Oxygenase into a pPET16B" and press ✨ enter ✨
"""


@click.command()
def main():
    print(WELCOME)
    while True:
        instruction = input(">")
        if instruction == "exit" or instruction == "quit" or instruction == "q":
            print("Goodbye 👋")
            break
        else:
            pass


if __name__ == "__main__":
    main()
