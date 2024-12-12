function data = jsoncDecode(jsoncString)
    % Remove any comments from the given JSON-C string and pass it to MATLAB's jsondecode
    data = jsondecode(jsonc.removeJSONComments(jsoncString));
end
