from langchain.agents import AgentExecutor, ZeroShotAgent
from langchain.agents.openai_functions_agent.base import OpenAIFunctionsAgent
from langchain.chat_models import ChatOpenAI

from ..tools import make_tools
from .prompt import FORMAT_INSTRUCTIONS


class AgentType:
    valid_models = {
        "ReactAgent": ZeroShotAgent,
        "OpenAIFunctionsAgent": OpenAIFunctionsAgent,
    }

    @classmethod
    def get_agent(cls, model_name: str = "ReactAgent"):
        return cls.valid_models[model_name]


class SanCloneAgent:
    def __init__(
        self,
        tools=None,
        llm=None,
        openai_api_key=None,
        temp=0.1,
        agent_type: str = "OpenAIFunctionsAgent",
        verbose=True,
    ):
        llm = ChatOpenAI(temperature=0.0, model="gpt-4", client=None)

        tools = make_tools(llm)
        self.agent_instance = AgentExecutor.from_agent_and_tools(
            tools=tools,
            agent=AgentType.get_agent(agent_type).from_llm_and_tools(llm, tools),
            return_intermediate_steps=True,
            handle_parsing_errors=True,
        )

    def run(self, prompt: str):
        return self.agent_instance.run(FORMAT_INSTRUCTIONS.format(input=prompt))
