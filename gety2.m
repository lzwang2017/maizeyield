function [gm, y2] = gety2(DPG, y1, y, WA)
  [n, p] = size(DPG);
  x = repmat(y1,1,p).*DPG;

  GM_A = zeros(0, p);
  GM_b = [];
  %% old
  % GM_Aeq = zeros(14, p);
  % GM_Aeq(:, [16:23, 73:78]) = eye(14);
  % GM_Aeq(:, 78+[16:23, 73:78]) = -eye(14);
  % GM_beq = zeros(14, 1);
  %% old
  %% new
  GM_Aeq = zeros(0, p);
  GM_beq = [];

  % GM_Aeq = zeros(65, p);
  % GM_Aeq(:, [16:35, 127:171]) = eye(65);
  % GM_Aeq(:, 171+[16:35, 127:171]) = -eye(65);
  % GM_beq = zeros(65, 1);
  %% new

  gm = zeros(p, 1);
  % rs = WA'.^2 * (x*gm-y).^2 / sum(WA.^2);

  fnall = [];
  fpall = [];

  while 1
    gm0 = QRGWAb(x, y, WA, GM_A, GM_b, -3e3*ones(p,1), 3e3*ones(p,1), [], GM_Aeq, GM_beq);
    % r0 = WA.^2 * (x*gm0 - y).^2 / sum(WA.^2);

    if isempty(gm0)
      break
    else
      gm = gm0;
    end

    fn = setdiff(find(x*gm < -1), fnall); % avoid same constraint added repeatedly due to computtaional errors.
    fp = setdiff(find(x*gm > 301), fpall);
    nfn = length(fn)
    nfp = length(fp)
    fnall = [fnall; fn];
    fpall = [fpall; fp];

    if isempty([fn; fp])
      break
    else
      GM_A = [GM_A; -x(fn,:); x(fp,:)];
      GM_b = [GM_b; -1 + zeros(nfn,1); 300*ones(nfp,1)];
    end
  end

  y2 = DPG * gm;
  % if ~isempty(gm)
  %   y2 = DPG * gm;
  % end
end
