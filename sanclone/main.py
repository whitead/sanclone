# -*- coding: utf-8 -*-
import click

WELCOME = """
Welcome to San Clone 👋  a molecular cloning agent 🧬.
Give it an instruction like "Clone NADH Oxidase from Streptococcus pyogenes into pET16b"
and press ✨ enter ✨
"""


@click.command()
def main():
    # check openai key
    try:
        from langchain.llms import OpenAI

        OpenAI(model="babbage-002")
    except Exception as e:
        if "OPENAI_API_KEY" in str(e):
            print("You need to set your OPENAI_API_KEY environment variable.")
            print("You can get one from https://beta.openai.com/")
            print("Then run the following command:")
            print("export OPENAI_API_KEY=<your key>")
            print("You can add this to your ~/.bashrc or ~/.bash_profile")
            return
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
