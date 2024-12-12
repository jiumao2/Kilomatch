function encodedString = jsonEncode(data)
    % Return a string representing the given data in the JSON standard
    try
        encodedString = jsonencode(data, "PrettyPrint", true);
    catch err
        % In Matlab R2020b and older, the PrettyPrint option is not
        % supported

        if not(err.identifier == "MATLAB:json:UnmatchedParameter")
            rethrow(err);
        end

        % Encode without formatting and add some newlines for some
        % basic formatting
        encodedString = jsonencode(data);

        fprintf([...
            '\n[esps.utils.jsonEncode] Note: You are using a Matlab version'...
            'that does not fully support the jsonencode function.\n'...
            'Consider upgrading to Matlab R2021a or newer for better '...
            'formatting of the JSON output files.\n'...
        ]);
    end
end