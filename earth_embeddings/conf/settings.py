from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    project_name: str = "earth_embeddings"
    debug: bool = False
