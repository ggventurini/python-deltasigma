function dsdemo4(action)
% A GUI-based demonstration of delta-sigma with sound.

% Copyright (c) 2001 Richard Schreier (Analog Devices, Inc.)

%To Do:
% Change color of 'Go' button when a parameter is changed and 'Go' needs
% to be pressed.

if nargin<1,
    action='initialize';
end;

if strcmp(action,'initialize'),
    clc
    fprintf(1,'\t\t\tDelta-Sigma Audio Demo \n\n\n');
    fprintf(1,'Set source, modulator and decimation filter paramters, ');
    fprintf(1,'then click "Go."\n');
    fprintf(1,'Click on input or output waveforms to listen to them.\n');
    fprintf(1,'Click on the output spectrum to listen to the amplified error.\n\n\n');
    shh = get(0,'ShowHiddenHandles');
    set(0,'ShowHiddenHandles','on');
    dsdemo4fig;
    set(0,'ShowHiddenHandles',shh)
    return

elseif strcmp(action,'Go')
    set(gcf,'Pointer','watch');
	fig_h = get(gco,'Parent');
	children = get(fig_h,'Children');
	axis_h = findobj(children,'Type','axes');
	% Compute u and plot it.
	T = get(findobj(children,'Tag','Duration'),'UserData');
	Fs = get(findobj(children,'Tag','Fs'),'UserData');
	DecFact = get(findobj(children,'Tag','DecFact'),'UserData');
	SincOrder = get(findobj(children,'Tag','SincOrder'),'UserData');
	FsOut = Fs/DecFact;
	SourceType = get(findobj(children,'Tag','SourceType'),'Value');
	N = round(T*Fs);
	switch SourceType
	case 1 	% sine
		SineAmp = get(findobj(children,'Tag','SineAmp'),'UserData');
		SineFreq = get(findobj(children,'Tag','SineFreq'),'UserData');
		u = SineAmp*sin(2*pi*SineFreq/Fs*[0:N-1]).*ds_hann(N);
		u0 = u(1:DecFact:end);
	case 2  % ramp
		u = linspace(-0.7,0.7,N);
		u0 = u(1:DecFact:end);
	case 3  % spech
		load(mfilename);	% for sd, ds
		u0 = ds;
		u = interp(sd,DecFact);
		N = length(u);
	end
	axes(findobj(axis_h,'Tag','u(t)'));
	hold on; cla;
	t = [0:length(u0)-1]/FsOut;
	h = plot(t,u0);
	set(h,'HitTest','off');
	set(gca,'UserData',u0);
	% Plot U(f), from 0 to FsOut/2
	axes(findobj(axis_h,'Tag','U(f)'));
	hold on; cla;
	if SourceType == 1 | SourceType == 3 
		N = length(u0);
		if SourceType == 1
			U = fft(u0)/(N/4);
		else
			U = fft(u0.*ds_hann(N))/(N/4);
		end
		f = linspace(0,FsOut,N+1); f = f(1:N/2+1);
		h=plot(f,dbv(U(1:N/2+1)));
		set(h,'HitTest','off');
	end
	% Compute v
    mod_h = findobj(get(fig_h,'Children'),'Tag','ModulatorFrame');
	ABCD = get(mod_h,'UserData');
	[v junk1 junk2 y] = simulateDSM(u,ABCD);
	clear junk1 junk2;
	q = v - y;		% Quantization error; assumes quantizer gain = 1.
	N = length(v);
	nPlot = 100;
	if N>nPlot
		n = floor( N/2-nPlot/2:N/2+nPlot/2-1 );
	else
		n = 1:N;
	end
	axes(findobj(axis_h,'Tag','v(t)'));
	hold on; cla;
    h = stairs(0:length(n)-1,v(n));
	set(h,'HitTest','off');
    h = plot(0:length(n)-1,u(n),'g');
	set(h,'HitTest','off');
	% Plot V(f), from 0 to Fs/2. Use the middle Nfft points of v
	N = length(v);
	Nfft = min( [N 16*8192] );
	n = ((N-Nfft)/2+1):((N+Nfft)/2);
	V = fft(v(n).*ds_hann(Nfft))/(Nfft/4);
	if SourceType==1
		inBin = round(SineFreq/Fs*Nfft+1);
	else
		inBin = ceil(Nfft/1000);
	end
	
	[f Vp] = logsmooth(V,inBin);
	axes(findobj(axis_h,'Tag','V(f)'));
	hold on; cla;
	h = semilogx( f*Fs, Vp );
	set(h,'HitTest','off');
	set(gca,'XLim',[100 Fs/2] );
	msg = sprintf('NBW = %.1f Hz ', Fs*1.5/Nfft);
	text(Fs/2,-90,msg,'Hor','right','Vert','mid')
	clear u V;
	% Compute w
	w = sinc_decimate(v,SincOrder,DecFact); clear v;
	filtered_q = sinc_decimate(q,SincOrder,DecFact); clear q;
	N = length(w);
	axes(findobj(axis_h,'Tag','w(t)'));
	hold on; cla;
	t = [0:length(w)-1]/FsOut;
	h = stairs(t,w);
	set(h,'HitTest','off');
	set(gca,'UserData',w);
    set(gcf,'Pointer','arrow');
	% Plot W(f), from 0 to FsOut/2
	axes(findobj(axis_h,'Tag','W(f)'));
	hold on; cla;
	set(gca,'UserData',filtered_q);
	if SourceType == 1 | SourceType == 3 
		Nfft = length(w);
		if SourceType == 1
			W = fft(w)/(N/4);
		else
			W = fft(w.*ds_hann(N))/(N/4);
		end
		f = linspace(0,FsOut,Nfft+1); f = f(1:Nfft/2+1);
		h=plot(f,dbv(W(1:Nfft/2+1)));
		msg = sprintf('  NBW = %.1f Hz', FsOut*1.5/Nfft);
		text(10,-90,msg,'Hor','left','Vert','mid')
		set(h,'HitTest','off');
	end

elseif strncmp(action,'mod',3) % mod1 or mod2
	btnParent_h = get(gco,'Parent');
	btnGroup_h = findobj(get(btnParent_h,'Children'),'Style','radiobutton');
	btn_h = findobj( btnGroup_h, 'Tag', action );
	for h = btnGroup_h'
		if h == btn_h
			set(h,'Value',1);
		else	
			set(h,'Value',0);
		end
	end
    mod_h = findobj(get(btnParent_h,'Children'),'Tag','ModulatorFrame');
	switch( action )
	case 'mod1'
		set(mod_h, 'UserData', [1 1 -1; 1 0 0]);
	case 'mod2'
		set(mod_h, 'UserData', [1 0 1 -1; 1 1 1 -2; 0 1 0 0]);
	end

elseif strcmp(action,'PlayU')
	u0 = get( gco, 'UserData' );
	FsOut = 8192;	% Output sample rate is fixed
	sound(u0,FsOut);

elseif strcmp(action,'PlayW')
	w = get( gco, 'UserData' );
	FsOut = 8192;	% Output sample rate is fixed
	sound(w,FsOut);
	
elseif strcmp(action,'PlayE')
	w = get( gco, 'UserData' );
	FsOut = 8192;	% Output sample rate is fixed
	soundsc(w,FsOut);

elseif strcmp(action,'Duration')
	h = gco;
	s = get(h,'String');
    v = get(h,'Userdata');
    T = round(eval(s,num2str(v)));
	set(h,'Userdata',T);
    set(h,'String', num2str(T,'%.0f'));

elseif strcmp(action,'Fs')
	h = gco;
	s = get(h,'String');
    v = get(h,'Userdata');
    Fs = 1e3*eval(s,num2str(v/1e3));
	FsOut = 8192;	% Output sample rate is fixed
	% Round Fs to a multiple of FsOut
	df = round(Fs/FsOut);
	df = max([1 df]);
	Fs = FsOut*df;
    set(h,'Userdata', Fs);
    set(h,'String', num2str(Fs/1e3,'%.0f'));
	% Change Decimation Factor
	h = findobj(gcf,'Tag','DecFact');
	set(h,'Userdata',df);
    set(h,'String', num2str(df,'%.0f'));

elseif strcmp(action,'DecFact')
	h = gco;
	s = get(h,'String');
    v = get(h,'Userdata');
    df = round(eval(s,num2str(v)));
	set(h,'Userdata',df);
    set(h,'String', num2str(df,'%.0f'));
	% Change Fs
	h = findobj(gcf,'Tag','Fs');
	FsOut = 8192;	% Output sample rate is fixed
	Fs = FsOut*df;
    set(h,'Userdata', Fs);
    set(h,'String', num2str(Fs/1e3,'%.0f'));
	
elseif strcmp(action,'SincOrder')
	h = gco;
	s = get(h,'String');
    v = get(h,'Userdata');
    SincOrder = round(eval(s,num2str(v)));
	set(h,'Userdata',SincOrder);
    set(h,'String', num2str(SincOrder,'%2d'));

elseif strcmp(action,'SourceType')	% sine, ramp, speech
	v = get(gco,'Value');
	h_list = get( get(gco,'Parent'), 'Children' );
	switch v
		case 1		%sine: Make SineFreq and SineAmp fields visible 
			visibility = 'on';
		case {2, 3}	%ramp or speech: Make SineFreq and SineAmp fields visible
			visibility = 'off';
	end
	h = findobj(h_list,'Tag','FreqLabel');
	set(h, 'Visible', visibility);
	h = findobj(h_list,'Tag','SineFreq');
	set(h, 'Visible', visibility);
	h = findobj(h_list,'Tag','AmpLabel');
	set(h, 'Visible', visibility);
	h = findobj(h_list,'Tag','SineAmp');
	set(h, 'Visible', visibility);
	
elseif strcmp(action,'SineFreq')
	h = gco;
	s = get(h,'String');
    v = get(h,'Userdata');
    SineFreq = round(eval(s,num2str(v)));
 	FsOut = 8192;	% Output sample rate is fixed
	if SineFreq >= FsOut/2
        waitfor(msgbox({'Anything above FsOut/2 will be inaudible.'},...
                 'Moddemo Error','error','modal'))
        SineFreq = v;
    end
    set(h,'Userdata',SineFreq,'String',num2str(SineFreq));

elseif strcmp(action,'SineAmp')
	h = gco;
	s = get(h,'String');
    v = get(h,'Userdata');
    SineAmp = eval(s,num2str(v));
	SineAmp = min([SineAmp 1]);
	SineAmp = max([SineAmp 0]);
    set(h,'Userdata',SineAmp,'String',num2str(SineAmp));
	
end    % if strcmp(action, ...lose all
