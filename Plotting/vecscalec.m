function  vecscalec(x,y,Img, disprange,xmap,ymap,vxmap,vymap,scl)
%VECSHOW is a modification of IMSHOW3D which displays 3D grayscale images in slice by slice fashion with mouse
%based slice browsing and window and level adjustment control. VECSHOW is
%intented for 2D+t images with vector maps.
%
% Usage:
% vecshow ( Image , [LOW HIGH], xmap, ymap, vxmap, vymap, scl )
%   
%    Image:      3D image MxNxK (K slices of MxN images) 
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the widest available range)
%    xmap:       matrix (time step, grid value Y, grid value X) with x positions for quiver plot
%    ymap:       matrix (time step, grid value Y, grid value X) with y positions for quiver plot
%    vxmap:      matrix (time step, grid value Y, grid value X) with x velocities for quiver plot
%    vymap:      matrix (time step, grid value Y, grid value X) with y velocities for quiver plot
%    scl:        scaling factor for quiver plot (-1 = no scaling)
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment). 
% 
% "Auto W/L" button adjust the window and level automatically 
%
% While "Fine Tune" check box is checked the window/level adjustment gets
% 16 times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse based window and level adjustment is set
% based on the user defined display intensity range; the wider the range
% the more sensitivity to mouse drag.
% 
% 
%   Example
%   --------
%       % Display an image (MRI example)
%       load mri 
%       Image = squeeze(D); 
%       figure, 
%       imshow3D(Image) 
%
%       % Display the image, adjust the display range
%       figure,
%       imshow3D(Image,[20 100]);
%
%   See also IMSHOW.

%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% 

sno = size(Img,3);  % number of slices
S = round(sno/2);

global InitialCoord;

MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end    

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
elseif numel(disprange) == 0
    [Rmin Rmax] = WL2R(Win, LevV);
else
    LevV = (double(disprange(2)) + double(disprange(1))) / 2;
    Win = double(disprange(2)) - double(disprange(1));
    WLAdjCoe = (Win + 1)/1024;
    [Rmin Rmax] = WL2R(Win, LevV);
end

FigPos = get(gcf,'Position');
S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [50 20 60 20];
Wval_Pos = [110 20 60 20];
Ltxt_Pos = [175 20 45 20];
Lval_Pos = [220 20 60 20];
BtnStPnt = uint16(FigPos(3)-250)+1;
if BtnStPnt < 300
    BtnStPnt = 300;
end
Btn_Pos = [BtnStPnt 20 100 20];
ChBx_Pos = [BtnStPnt+110 20 100 20];

if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end    
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);

set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)

% ax = gca;
mag = sqrt(vxmap.^2 + vymap.^2);
magmax = prctile(mag(:),99);
vxmap = vxmap./magmax;
vymap = vymap./magmax;

dx = abs(xmap(1,2)-xmap(1,1));
dy = abs(ymap(2,1)-ymap(1,1));
dgrid = sqrt(dx^2+dy^2);
vxmap = vxmap.*(dgrid*2); % longest vector is 2 times the grid spacing
vymap = vymap.*(dgrid*2); % longest vector is 2 times the grid spacing
        
hold on;
% quivercEdit(xmap(:,:,S),ymap(:,:,S),vxmap(:,:,S),vymap(:,:,S),scl, umax,ax);
quivax = axes;
% set(quivax, 'Position', imax.Position)
% quivax.Visible = 'off';

nonNaNind = find(~isnan(vr(:)));
xyuvvrNN = [x(nonNaNind),y(nonNaNind),u(nonNaNind),v(nonNaNind),vr(nonNaNind)];
[xyuvvrNNs, ~] = sortrows(xyuvvrNN,5);
iCCs=round(linspace(1,size(CC,1),n));
iData = round(linspace(0,size(xyuvvrNNs,1),n+1));

xs = x(1);
xf = x(end);
xa = (3 * xs + xf)/4;
xb = (3 * xf + xs)/4;
ys = y(1);
yf = y(end);
ya = (3 * ys + yf)/4;
yb = (3 * yf + ys)/4;
axes(ax)
h=imagesc([xa,xb],[ya,yb],reshape(xyuvvrNNs([1,end;1,end],5),2,2));
set(h,'visible','off');
caxis(xyuvvrNNs([1,end],5));
hold on;
for it=1:n
    c = CC(iCCs(it),:);
    si = iData(it)+1;
    ei = iData(it+1);
    hh=quiver(xyuvvrNNs(si:ei,1),xyuvvrNNs(si:ei,2),xyuvvrNNs(si:ei,3),xyuvvrNNs(si:ei,4),scale*it/n,'Color',c);
end    
colorbar(ax);

q = quiverc2wcmap(xmap,ymap, vxmap(:,:,S),vymap(:,:,S), scl, quivax);


imax = axes('position',[0.1,0.25,0.8,0.7]), imagesc(x,y,Img(:,:,S), [Rmin Rmax]);
colormap gray
axis image

hold off;


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-250)+1;
        if BtnStPnt < 300
            BtnStPnt = 300;
        end
        Btn_Pos = [BtnStPnt 20 100 20];
        ChBx_Pos = [BtnStPnt+110 20 100 20];
        if sno > 1
            set(shand,'Position', S_Pos);
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        h = get(gca,'children');
        delete(findobj(gca,'Type','Quiver'));
        set(findobj(h,'Type','Image'),'cdata',Img(:,:,S))
      
        hold on;
        % quivercEdit(xmap(:,:,S),ymap(:,:,S),vxmap(:,:,S),vymap(:,:,S),scl, umax,ax);
        quiver(xmap,ymap,vxmap(:,:,S),vymap(:,:,S),scl,'g');
        hold off;
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        h = get(gca,'children');
        delete(findobj(gca,'Type','Quiver'));
        
        set(findobj(h,'Type','Image'),'cdata',Img(:,:,S))
        
        hold on;
        % quivercEdit(xmap(:,:,S),ymap(:,:,S),vxmap(:,:,S),vymap(:,:,S),scl, umax,ax);
        quiver(xmap,ymap,vxmap(:,:,S),vymap(:,:,S),scl,'g');
        hold off;
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

end
% -=< Maysam Shahedi (mshahedi@gmail.com), April 19, 2013>=-