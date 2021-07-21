function [gm, ap, yh] = train_Dec9(tpower, uclusters, vmem, RT, DPG, WA, y, gm, ap)
fig01 = 1;
tic
ptr = zeros(100, 2);
ip = 0;
n = length(vmem);
nmem = length(unique(vmem));
nfT = round(0.8 * n);
nfV = n - nfT;

if ~exist('gm', 'var') || ~exist('ap', 'var')
    y1 = y/max(y);
    [gm, y2] = gety2(DPG, y1, y, WA);
    [ap, y1] = gety1(tpower, uclusters, vmem, RT, y2, y, WA);
else
    y1 = RT * ap;
    y2 = DPG * gm;
end

rs = sqrt(sum(WA.^2.*(y-y1.*y2).^2)/sum(WA.^2))

ok0 = 2;
ok = ok0;
ip = ip + 1;
ptr(ip, :) = [toc, rs];
while ok >= 1
    while 1
        fT = randperm(n);
        fV = fT(nfT+[1:nfV]);
        fT = fT(1:nfT);
        if length(unique(vmem(fT))) == nmem
            break
        end
    end

    gm0 = gety2(DPG(fT, :), y1(fT), y(fT), WA(fT));

    if sum(abs(gm0)) < 1e-4
        gm0 = gm;
        rT = sqrt(sum(WA(fT).^2.*(y(fT)-y1(fT).*y2(fT)).^2)/sum(WA(fT).^2));
        rV = sqrt(sum(WA(fV).^2.*(y(fV)-y1(fV).*y2(fV)).^2)/sum(WA(fV).^2));
        rT0 = rT;
        rV0 = rV;
        r2 = rs;
    else
        gm0 = 0.8*gm + 0.2*gm0;
        y02 = DPG * gm0;
        rT = sqrt(sum(WA(fT).^2.*(y(fT)-y1(fT).*y2(fT)).^2)/sum(WA(fT).^2));
        rV = sqrt(sum(WA(fV).^2.*(y(fV)-y1(fV).*y2(fV)).^2)/sum(WA(fV).^2));
        rT0 = sqrt(sum(WA(fT).^2.*(y(fT)-y1(fT).*y02(fT)).^2)/sum(WA(fT).^2));
        rV0 = sqrt(sum(WA(fV).^2.*(y(fV)-y1(fV).*y02(fV)).^2)/sum(WA(fV).^2));
        r2 = sqrt(sum(WA.^2.*(y-y1.*y02).^2)/sum(WA.^2));
    end

    [rT0 rV0 r2; rT rV rs]
    % [r2 rs]
    if rT0 <= rT && rV0 <= rV && r2 <= rs && (rT0+rV0+r2) < (rT+rV+rs) - 1e-3
        % if r2 <= rs
        rs = r2;
        gm = gm0;
        y2 = y02;
        ok = ok0;
    else
        ok = ok - 1;
    end
    ip = ip + 1;
    ptr(ip, :) = [toc, r2];
    if fig01 == 1
        plot(ptr(1:ip,1), ptr(1:ip,2));
        drawnow
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while 1
        fT = randperm(n);
        fV = fT(nfT+[1:nfV]);
        fT = fT(1:nfT);
        if length(unique(vmem(fT))) == nmem
            break
        end
    end

    ap0 = gety1(tpower, uclusters, vmem(fT), RT(fT, :), y2(fT), y(fT), WA(fT));

    ap0 = 0.8*ap + 0.2*ap0;
    y01 = RT * ap0;
    rT = sqrt(sum(WA(fT).^2.*(y(fT)-y1(fT).*y2(fT)).^2)/sum(WA(fT).^2));
    rV = sqrt(sum(WA(fV).^2.*(y(fV)-y1(fV).*y2(fV)).^2)/sum(WA(fV).^2));
    rT0 = sqrt(sum(WA(fT).^2.*(y(fT)-y01(fT).*y2(fT)).^2)/sum(WA(fT).^2));
    rV0 = sqrt(sum(WA(fV).^2.*(y(fV)-y01(fV).*y2(fV)).^2)/sum(WA(fV).^2));
    r1 = sqrt(sum(WA.^2.*(y-y01.*y2).^2)/sum(WA.^2));
    % [r1 rs]
    [rT0 rV0 r1; rT rV rs]
    if rT0 <= rT && rV0 <= rV && r1 <= rs && (rT0+rV0+r1) < (rT+rV+rs) - 1e-3
        % if r1 <= rs
        rs = r1;
        ap = ap0;
        y1 = y01;
        ok = ok0;
    else
        ok = ok - 1;
    end
    ip = ip + 1;
    ptr(ip, :) = [toc, r1];
    if fig01 == 1
        plot(ptr(1:ip,1), ptr(1:ip,2));
        drawnow
    end
end

ptr(1:ip,:)
yh = y1 .* y2;
end
