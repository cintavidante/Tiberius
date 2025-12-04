from astropy.io import fits

def get_fits_image_extensions_from_config(instr_cfg, filename):
    rule_block = instr_cfg.image_extension
    rule = rule_block["rule"]

    if rule == "explicit":
        return rule_block["list"]
    
    # This logic is for ACAM data: HDU0 = primary, windows start at HDU1
    elif rule == "nwin_map":
        with fits.open(filename) as hdul:
            n_hdus = len(hdul)

        nwin = n_hdus - 1

        if nwin not in rule_block["map"]:
            raise ValueError(f"nwin={nwin} not in the map for instruments={instr_cfg.name}")
        
        return rule_block["map"][nwin]
    
    else:
        raise ValueError(f"Unknown rule {rule}")
    

def compute_buffer_pixels(instr_cfg, row, oversampling_factor):
    """Compute buffer pixels for left and right edges based on YAML-defined rules"""

    rule = instr_cfg.get("buffer_rule", {"type": "static"})

    # Static rule (ACAM, EFOSC, JWST, etc.)
    if rule["type"] == "static":
        base_left = rule.get("left", 0)
        base_right = rule.get("right", 0)
        return base_left * oversampling_factor, base_right * oversampling_factor

    # KECK rule
    if rule["type"] == "edge_threshold":
        left_thresh = rule["left_threshold"]
        right_thresh = rule["right_threshold"]
        buf = rule["buffer_pixels"]
        default = rule["default_buffer_pixels"]
        mul = rule.get("multiply_by_oversampling", False)

        # Compute dynamic buffers
        left = buf if row[0] > left_thresh else default
        right = buf if row[-1] > right_thresh else default

        if mul:
            left *= oversampling_factor
            right *= oversampling_factor

        return left, right

    raise ValueError(f"Unknown buffer_rule type: {rule['type']}")
