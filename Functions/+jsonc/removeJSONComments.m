function parsedString = removeJSONComments(subjectString)
    % Removes comments from a JSON-C string and returns a valid JSON string
    % Valid comments are:
    %  - One-line comments: // ....
    %  - Multi-line comments: /* ... */
    %  - Comment delimiters are ignored if they are inside "" strings
    % This function uses a regular expression to remove the comments,
    %   see documentation on MATLAB-style regular expressions.

    % This regexp matches strings delimited by " as well as comments
    % This is to ensure that we don't match // and /**/ delimiters inside strings
    pattern = '("[^"]*?")|(/\*.+?\*/)|(//.+?$)';
    [match,noMatch] = regexp(subjectString, pattern, 'match', 'split', 'lineanchors');

    % Iterate though the matches and only remove the ones that are not strings
    for k = 1:length(match)
        if not(startsWith(match{k}, '"'))
            match{k} = '';
        end
    end

    % Reassemble the result
    parsedString = strjoin(noMatch, match);
end
