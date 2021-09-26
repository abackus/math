% Modify mod to avoid off-by-one bugs
% Intead of using 0 as a rep, uses the order k of the cyclic group as a rep

function l = obob_mod(n, k)
    l = mod(n, k);
    if l == 0
        l = k;
    end
end
