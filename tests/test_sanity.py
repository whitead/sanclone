import sanclone


def test_version():
    assert sanclone.__version__

def test_echo_tool():
    from sanclone.tools import EchoTool
    from sanclone import State
    tool = EchoTool(shared_state=State())
    assert tool.run("Hello") == "Hello"
    