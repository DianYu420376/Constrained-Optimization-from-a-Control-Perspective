function arr = insertElement(arr, a, p)
    arr = [arr(1:p-1); a; arr(p:end)];
end