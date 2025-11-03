def print_dict_mismatches(d1, d2, rtol=1e-6, atol=1e-8):
    """compare two dicts and print mismatched fields with details.
    """
    k1, k2 = set(d1.keys()), set(d2.keys())
    mismatches = []

    # missing keys
    for k in sorted(k1 - k2):
        mismatches.append((k, d1.get(k), '<missing>'))
    for k in sorted(k2 - k1):
        mismatches.append((k, '<missing>', d2.get(k)))

    import numpy as np

    def is_array_like(x):
        return isinstance(x, (list, tuple, np.ndarray))

    def to_numeric_array(x):
        """Try convert to float array. Return (array, is_numeric)."""
        try:
            a = np.asarray(x, dtype=float)
            return a, True
        except Exception:
            return np.asarray(x, dtype=object), False

    # common keys with different values
    for k in sorted(k1 & k2):
        v1, v2 = d1[k], d2[k]
        equal = False
        if is_array_like(v1) or is_array_like(v2):
            a1, n1 = to_numeric_array(v1)
            a2, n2 = to_numeric_array(v2)
            if a1.shape != a2.shape:
                equal = False
            else:
                if n1 and n2:
                    equal = np.allclose(a1, a2, rtol=rtol, atol=atol, equal_nan=True)
                else:
                    # 非数值数组：用元素相等比较
                    eq = np.equal(a1, a2)
                    eq = np.asarray(eq)
                    equal = bool(np.all(eq))
        elif isinstance(v1, float) or isinstance(v2, float):
            try:
                equal = np.isclose(float(v1), float(v2), rtol=rtol, atol=atol, equal_nan=True)
            except Exception:
                equal = (v1 == v2)
        else:
            equal = (v1 == v2)
        if not equal:
            mismatches.append((k, v1, v2))

    def summarize(v):
        if isinstance(v, np.ndarray):
            if v.size == 0:
                return 'array([], shape=(0,))'
            if np.issubdtype(v.dtype, np.number):
                return f'array(shape={v.shape}, min={np.nanmin(v):.6g}, max={np.nanmax(v):.6g})'
            # 非数值数组
            flat = v.ravel()
            head = flat[:5].tolist()
            tail = flat[-3:].tolist() if flat.size > 8 else []
            if tail:
                return f'array(dtype={v.dtype}, shape={v.shape}, head={head} ... tail={tail})'
            return f'array(dtype={v.dtype}, shape={v.shape}, data={flat.tolist()})'
        if isinstance(v, (list, tuple)):
            lv = list(v)
            if len(lv) > 8:
                return f'{type(v).__name__}(len={len(lv)}, head={lv[:5]} ... tail={lv[-3:]})'
            return f'{type(v).__name__}({lv})'
        if isinstance(v, float):
            return f'{v:.10g}'
        return repr(v)

    print(f'Mismatched fields ({len(mismatches)}):')
    for k, v1, v2 in mismatches:
        if is_array_like(v1) or is_array_like(v2):
            a1, n1 = to_numeric_array(v1)
            a2, n2 = to_numeric_array(v2)
            if a1.shape == a2.shape and a1.size > 0:
                if n1 and n2:
                    diff = np.abs(a1 - a2)
                    notclose = ~np.isclose(a1, a2, rtol=rtol, atol=atol, equal_nan=True)
                    cnt = int(np.count_nonzero(notclose))
                    mad = float(np.nanmax(diff)) if diff.size else float('nan')
                    print(f'  - {k}: numeric array, shape={a1.shape}, not_close={cnt}, max_abs_diff={mad:.6g}')
                else:
                    eq = np.equal(a1, a2)
                    cnt = int(a1.size - int(np.count_nonzero(eq)))
                    # 打印部分不等样本
                    idx = np.where(~eq.ravel())[0]
                    samples = []
                    for i in idx[:5]:
                        samples.append((a1.ravel()[i], a2.ravel()[i]))
                    extra = ' ...' if cnt > 5 else ''
                    print(f'  - {k}: non-numeric array, shape={a1.shape}, unequal={cnt}, samples={samples}{extra}')
            else:
                print(f'  - {k}: shape_ref={getattr(np.asarray(v1, dtype=object), "shape", None)} != shape_wr={getattr(np.asarray(v2, dtype=object), "shape", None)}')
        else:
            print(f'  - {k}: ref={summarize(v1)}, wr={summarize(v2)}')
    return mismatches