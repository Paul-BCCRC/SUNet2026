function CB2EXCEL_function( inpath, infilename, ext )

    % [ features ] = readCBX(inpath , infilename, ext );
    % [~, F] = size(features);
    % filename = [infilename,'_',ext,'.xlsx'];
    % for f=1:F
    %     writecell(features(1, f), filename, "Sheet", 1, "Range", ['A',num2str(f)]);
    % end
    % writecell(features(2:end, :)', filename, "Sheet", 1, "Range", 'B1');
    % movefile(filename, inpath);

    [ features ] = readCBX(inpath , infilename, ext );
    [~, F] = size(features);
    filename = [infilename,'_',ext,'.xlsx'];
    for f=1:F
        writecell(features(1, f), filename, "Sheet", 1, "Range", [num2letters(f),'1']);
    end
    writecell(features(2:end, :), filename, "Sheet", 1, "Range", 'A2', 'AutoFitWidth', false); 

    % hExcel = actxserver('Excel.Application');
    % hWorkbook = hExcel.Workbooks.Open(fullfile(pwd, filename));
    % hWorksheet = hWorkbook.Sheets.Item(1);
    % 
    % % Set first column width to 30
    % hWorksheet.Columns.Item(1).ColumnWidth = 30; 
    % 
    % hWorkbook.Save;
    % hExcel.Quit;


    movefile(filename, inpath);
end