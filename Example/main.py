from metabomix import MetaboMix

def main():
    settings_path: str = "/lustre/BIF/nobackup/hendr218/metabomix/Example/Example_recipe.json"
    example_mix = MetaboMix(settings_path)
    example_mix.run_all()
if __name__ == "__main__":
    main()