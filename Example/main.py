from metabomix import MetaboMix

def main() -> None:
    settings_path: str = "recipe.json"
    example_mix = MetaboMix(settings_path)
    example_mix.run_all()

if __name__ == "__main__":
    main()