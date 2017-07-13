if(isnan(fmax)&& isnan(fmin))
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'Padding',0); 
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'Padding',0); 
                    end            
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc,'Padding',0); 
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'f0',fc,'Padding',0); 
                    end   
            end
        elseif(isnan(fmax))
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'Padding',0); 
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'Padding',0); 
                    end
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc,'Padding',0); 
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'f0',fc,'Padding',0); 
                    end
            end
        elseif(isnan(fmin))
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'Padding',0); 
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'Padding',0); 
                    end
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc,'Padding',0); 
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'f0',fc,'Padding',0); 
                    end
            end
        else
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'Padding',0);
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'Padding',0); 
                    end
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc,'Padding',0);
                    else
                        [WT,handles.freqarr]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Window',wavelet_type_selected,'f0',fc,'Padding',0); 
                    end                 
            end
        end