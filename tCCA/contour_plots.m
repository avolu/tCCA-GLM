function [] = contour_plots(METRIC, ttl, evparams, pOpt, cntno, cmap, flip, lopttext)
%CONTOUR_PLOTS Summary of this function goes here
%   Detailed explanation goes here
for tl = 1:numel(evparams.stpsize)
xtl{tl} = num2str(evparams.stpsize(tl)/25, '%.2g');
xtl{tl}=strrep(xtl{tl}, '0.', '.');
end

%% create combined surface plots (depict objective function)
[X,Y] = meshgrid(evparams.stpsize,evparams.tlags);

%% Plot Objective function results
figure
climits = [min(METRIC(:)) max(METRIC(:))];
for ii=1:10
    subplot(2,5,ii)
    contourf(X,Y, squeeze(METRIC(:,:,ii)), cntno)
    xlabel('stepsize / s')
    xticks(evparams.stpsize(1:2:end))
    xticklabels(xtl(1:2:end))
    ylabel('time lags / s')
    title([ttl ', ct = ' num2str(evparams.cthresh(ii))])
    grid on
    buf =  squeeze(METRIC(:,:,ii));
    switch flip
        case 'max'
            colormap(cmap)
            [r,c] = ind2sub(size(buf),find(buf == max(buf(:))));
            limit = climits(2);
        case 'min'
            colormap(flipud(cmap))
            [r,c] = ind2sub(size(buf),find(buf == min(buf(:))));
            limit = climits(1);
    end
    if ~mod(ii,5)
        cb = colorbar;
        if ii/5 == 1
            cb.Position = [0.917027047511808,0.582968760552204,0.012213203762342,0.342857142957143]; %cb.Position + 1e-10;
        else
            cb.Position = [0.917027047511808,0.106986227801112,0.012213203762342,0.342857142957143]
        end
    end
    caxis(climits)
    % mark local optima
    hold on
    if squeeze(METRIC(r(1),c(1),ii)) == limit
            mrkc = 'g';
        else
            mrkc = 'k';
        end
        if numel(r)<4
            p=numel(r);
        else
            p=1;
        end
        for pp = 1:p
            plot(evparams.stpsize(c(pp)),evparams.tlags(r(pp)),'ko','MarkerFaceColor', mrkc)
            if pp==1
                if lopttext
                text(evparams.stpsize(c(pp)),evparams.tlags(r(pp)), ['\leftarrow ' num2str(METRIC(r(pp),c(pp),ii))])
                end
            end
        end
        % mark optimum from objective function
        if ii == pOpt(1,3)
            plot(evparams.stpsize(pOpt(1,2)),evparams.tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
            text(evparams.stpsize(pOpt(1,2)),evparams.tlags(pOpt(1,1)), ['\leftarrow ' num2str(METRIC(pOpt(1,1),pOpt(1,2),ii))])
        end
end
end

