from kml_rnaseq import get_threads_dict, get_conda_env_dict, get_database_dict


def test_threads():
    print(get_threads_dict())


def test_conda_env():
    print(get_conda_env_dict())


def test_database():
    print(get_database_dict())
