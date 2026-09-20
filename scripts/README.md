# Scripts

## Package map

Regenera el mapa navegable de módulos y símbolos públicos de `euskera/` con:

```bash
python scripts/generate_package_map.py
```

El script analiza los módulos con `ast`, omite imports y nombres privados y
respeta `__all__` cuando está definido. Ejecuta la comprobación sin escribir
ningún archivo con:

```bash
python scripts/generate_package_map.py --check
```

Después de añadir o quitar módulos o símbolos públicos, ejecuta el comando
normal y revisa `docs/package-map.md` junto con el resto de la documentación.
